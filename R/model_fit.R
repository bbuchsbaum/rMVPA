#' @keywords internal
#' @noRd
requireNamespaceQuietStop <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop(paste('package',package,'is required'), call. = FALSE)
}





#' @keywords internal
#' @noRd
get_control <- function(y, nreps) { # nreps for bootstrap tuning
  is_class <- is.factor(y)
  
  # Determine primary metric for tuning and optimization
  metric_name <- if (is_class) {
    if (nlevels(y) == 2) "roc_auc" else "accuracy" # Defaulting to accuracy for multiclass
  } else {
    "rmse" 
  }
  
  list(
    metric = metric_name, 
    number = nreps # For bootstrap reps in tuning, if applicable
  )
}


#'
#' This function finds the best hyperparameters for a given model specification
#' using a specified tuning grid and cross-validation.
#'
#' @param mspec A model specification derived from the \code{mvpa_model} class.
#' @param x The training data matrix.
#' @param y The response vector.
#' @param wts Optional class weights (if the underlying model supports it).
#' @param param A \code{data.frame} representing the tuning grid, where
#'        parameter names are indicated by column names.
#' @param nreps The number of bootstrap replications (default is 10).
#' @return A data frame containing the best hyperparameter values.
#' @keywords internal
#' @noRd
#' @importFrom rsample bootstraps
#' @importFrom purrr map_dfr map_dbl
#' @importFrom tibble tibble
tune_model <- function(mspec, x, y, wts, param_grid, nreps = 10) {
  
  # Get metric to optimize from the mspec's control object (derived via get_control)
  # The 'model' element within mspec is the list from MVPAModels
  control_obj <- get_control(y, nreps) 
  metric_to_optimize <- control_obj$metric 
  
  # Define the yardstick metric function based on metric_to_optimize
  # And define if higher is better for this metric
  higher_is_better <- TRUE
  metric_fn <- switch(metric_to_optimize,
    "roc_auc"     = yardstick::roc_auc_vec,
    "mn_log_loss" = { higher_is_better <- FALSE; yardstick::mn_log_loss_vec },
    "accuracy"    = yardstick::accuracy_vec,
    "rmse"        = { higher_is_better <- FALSE; yardstick::rmse_vec },
    "rsq"         = yardstick::rsq_vec,
    stop("Unsupported metric for tuning: ", metric_to_optimize)
  )
  
  y_vector <- if(is.matrix(y) && ncol(y) == 1) y[,1] else y
  
  df_for_rsample <- as.data.frame(x)
  df_for_rsample$.response_var_for_stratification <- y_vector
  
  # Using bootstraps for tuning
  resamples_obj <- if(is.factor(y_vector)) {
    rsample::bootstraps(df_for_rsample, times = nreps, strata = tidyselect::all_of(".response_var_for_stratification"))
  } else {
    rsample::bootstraps(df_for_rsample, times = nreps)
  }

  tuning_metrics <- purrr::map_dfr(seq_len(nrow(param_grid)), .f = function(param_idx) {
    current_params_df <- param_grid[param_idx, , drop = FALSE]
    
    # Performance over resamples for this parameter set
    resample_perf <- purrr::map_dbl(resamples_obj$splits, .f = function(split) {
      train_df_fold <- rsample::analysis(split)
      test_df_fold  <- rsample::assessment(split)
      
      y_train_fold <- train_df_fold$.response_var_for_stratification
      x_train_fold <- as.matrix(train_df_fold[, !(names(train_df_fold) %in% ".response_var_for_stratification"), drop = FALSE])
      
      y_test_fold  <- test_df_fold$.response_var_for_stratification
      x_test_fold  <- as.matrix(test_df_fold[, !(names(test_df_fold) %in% ".response_var_for_stratification"), drop = FALSE])

      # mspec$model is the list from MVPAModels (e.g., MVPAModels$sda_notune)
      # Call its $fit element
      # The `param` argument to the model's fit function should be the current set
      fit_obj <- mspec$model$fit(x_train_fold, y_train_fold, 
                                 wts = wts, # Pass weights if available
                                 param = current_params_df, 
                                 lev = levels(y_vector), 
                                 last = FALSE,  # Not the last model
                                 weights = NULL,  # Different from wts
                                 classProbs = is.factor(y_vector)) # Pass classProbs

      # Predict and evaluate
      metric_value <- NA_real_
      if (metric_to_optimize == "roc_auc" || metric_to_optimize == "mn_log_loss") {
         preds_probs <- mspec$model$prob(fit_obj, x_test_fold)
         if (is.factor(y_test_fold) && is.matrix(preds_probs) && !is.null(levels(y_test_fold))) {
            colnames(preds_probs) <- levels(y_test_fold)
         }
         
         # For binary classification with roc_auc, pass only the probability of the second class
         if (metric_to_optimize == "roc_auc" && nlevels(y_test_fold) == 2) {
           prob_positive <- preds_probs[, levels(y_test_fold)[2]]
           metric_value <- metric_fn(truth = y_test_fold, estimate = prob_positive, event_level = "second")
         } else {
           metric_value <- metric_fn(truth = y_test_fold, estimate = preds_probs)
         }
      } else {
         preds <- mspec$model$predict(fit_obj, x_test_fold)
         if (is.factor(y_test_fold) && !is.null(levels(y_test_fold))) {
            preds <- factor(preds, levels = levels(y_test_fold))
         }
         metric_value <- metric_fn(truth = y_test_fold, estimate = preds)
      }
      metric_value
    }) # End loop over resamples_obj$splits
    
    # Return mean performance for this parameter set
    data.frame(.param_id = param_idx, mean_metric = mean(resample_perf, na.rm = TRUE))
  }) # End loop over param_grid
  
  best_row_idx <- if (higher_is_better) {
    which.max(tuning_metrics$mean_metric)
  } else {
    which.min(tuning_metrics$mean_metric)
  }
  
  if (length(best_row_idx) == 0 || is.na(tuning_metrics$mean_metric[best_row_idx])) { 
    warning("tune_model: All parameter combinations resulted in NA performance. Returning first parameter set.")
    best_row_idx <- 1
  }
  
  best_param_id <- tuning_metrics$.param_id[best_row_idx]
  best_tune_params <- param_grid[best_param_id, , drop = FALSE]
  
  return(best_tune_params)
}

#' Fit an MVPA model
#'
#' This function fits a multivariate pattern analysis (MVPA) model to the given data.
#'
#' @param obj An object derived from the \code{mvpa_model} class.
#' @param x The training data matrix.
#' @param y The response vector.
#' @param wts Optional class weights (if the underlying model supports it).
#' @param param The hyperparameters of the model.
#' @param classProbs Logical; if TRUE, class probabilities should be computed (default is FALSE).
#' @param ... Additional arguments to be passed to the underlying model fitting function.
#' @return A fitted model object with additional attributes "obsLevels" and "problemType".
#' @rdname fit_model-methods
#' @noRd
fit_model.mvpa_model <- function(obj, x, y, wts, param, classProbs, ...) {
  fit <- obj$model$fit(x,y,wts=wts,param=param,lev=levels(y), classProbs=classProbs, ...)
  
  # Add levels both as attribute and list element for consistency
  fit$obsLevels <- levels(y)
  attr(fit, "obsLevels") <- levels(y)
  
  if (is.factor(y)) {
    attr(fit, "problemType") <- "Classification"
  } else {
    attr(fit, "problemType") <- "Regression"
  }
  
  fit
}

#' Predict class labels and probabilities for new data using a fitted model
#'
#' @param object A fitted model object of class \code{class_model_fit}.
#' @param newdata New data to predict on, either as a \code{matrix} or a \code{NeuroVec} or \code{NeuroSurfaceVector} object.
#' @param sub_indices The subset of row indices to compute predictions on (optional).
#' @param ... Additional arguments to be passed to the underlying prediction function.
#' @return A list containing class predictions and probabilities with class attributes "classification_prediction", "prediction", and "list".
#' @noRd
#' @keywords internal
predict.class_model_fit <- function(object, newdata, sub_indices=NULL,...) {
  tryCatch({
    mat <- if (inherits(newdata, "NeuroVec") || inherits(newdata, "NeuroSurfaceVector")) {
      series(newdata, object$fit$vox_ind)
    } else {
      newdata
    }
    
    if (!is.null(sub_indices)) {
      assert_that(is.vector(sub_indices))
      mat <- mat[sub_indices,,drop=FALSE]
    }
    
    if (!is.null(object$feature_mask)) {
      mat <- mat[, object$feature_mask,drop=FALSE]
    }

    futile.logger::flog.debug("Predicting with data dimensions: %s", paste(dim(mat), collapse=" x "))
    
    probs <- object$model$prob(object$fit, mat)
    if (is.null(probs) || length(probs) == 0) {
      stop("Model probability calculation returned NULL or empty result")
    }
    
    colnames(probs) <- levels(object$y)
    cpred <- max.col(probs)
    cpred <- levels(object$y)[cpred]
    ret <- list(class=cpred, probs=probs)
    class(ret) <- c("classification_prediction", "prediction", "list")
    ret
    
  }, error = function(e) {
    futile.logger::flog.error("Class model prediction failed: %s", e$message)
    futile.logger::flog.debug("Input data dimensions: %s", paste(dim(newdata), collapse=" x "))
    stop(sprintf("Prediction failed: %s", e$message))
  })
}



#' Predict continuous values for a new dataset using a regression model
#'
#' This function predicts continuous values for new data using a fitted regression model.
#'
#' @param object A fitted model object of class \code{regression_model_fit}.
#' @param newdata New data to predict on, either as a matrix or a \code{NeuroVec} or \code{NeuroSurfaceVector} object.
#' @param sub_indices A vector of indices used to subset rows of `newdata` (optional).
#' @param ... Additional arguments to be passed to the underlying prediction function.
#' @return A list containing predicted continuous values with class attributes "regression_prediction", "prediction", and "list".
#' @noRd
#' @keywords internal
predict.regression_model_fit <- function(object, newdata, sub_indices=NULL,...) {
  #browser()
  tryCatch({
    mat <- if (inherits(newdata, "NeuroVec") || inherits(newdata, "NeuroSurfaceVector")) {
      series(newdata, object$fit$vox_ind)
    } else {
      newdata
    }
    
    if (!is.null(sub_indices)) {
      assert_that(is.vector(sub_indices))
      mat <- mat[sub_indices,,drop=FALSE]
    }
    
    if (!is.null(object$feature_mask)) {
      mat <- mat[, object$feature_mask,drop=FALSE]
    }

    futile.logger::flog.debug("Regression prediction with data dimensions: %s", paste(dim(mat), collapse=" x "))
    
    preds <- object$model$predict(object$fit, mat)
    if (is.null(preds) || length(preds) == 0) {
      stop("Model prediction returned NULL or empty result")
    }
    
    ret <- list(preds=preds)
    class(ret) <- c("regression_prediction", "prediction", "list")
    ret
    
  }, error = function(e) {
    futile.logger::flog.error("Regression model prediction failed: %s", e$message)
    futile.logger::flog.debug("Input data dimensions: %s", paste(dim(newdata), collapse=" x "))
    stop(sprintf("Prediction failed: %s", e$message))
  })
}


#' @rdname merge_predictions-methods
#' @export
#' @method merge_predictions regression_prediction
merge_predictions.regression_prediction <- function(obj1, rest, ...) {
  args <- list(...)
  weights <- if (!is.null(args$weights)) args$weights else rep(1,length(rest)+1)/(length(rest)+1)

  allobj <- c(obj1, rest)
  assert_that(all(sapply(allobj, function(obj) inherits(obj, "regression_prediction"))))
  
  preds <- lapply(1:length(allobj), function(i) {
    allobj[[i]]$preds * weights[i]
  })
  
  final_pred <- rowMeans(do.call(cbind, preds))
  ret <- list(preds=final_pred)
  class(ret) <- c("regression_prediction", "prediction", "list")
  ret
}


#' @rdname merge_predictions-methods
#' @export
#' @method merge_predictions classification_prediction
merge_predictions.classification_prediction <- function(obj1, rest, ...) {
  args <- list(...)
  weights <- if (!is.null(args$weights)) args$weights else rep(1,length(rest)+1)/(length(rest)+1)

  allobj <- vector(mode="list", length(rest)+1)
  allobj[[1]] <- obj1
  allobj[2:length(allobj)] <- rest
  
  assert_that(all(sapply(allobj, function(obj) inherits(obj, "classification_prediction"))))
  
  preds <- lapply(1:length(allobj), function(i) {
    allobj[[i]]$prob * weights[i]
  })
  
  prob <- preds[!sapply(preds, function(x) is.null(x))]
  pfinal <- Reduce("+", prob)
  
  cnames <- colnames(pfinal)
  maxids <- apply(pfinal, 1, which.max)
  len <- sapply(maxids, length)
  
  if (any(len == 0)) {
    maxids[len == 0] <- NA
  }
  
  maxids <- unlist(maxids)
  
  pclass <- cnames[maxids]
  ret <- list(class=pclass, probs=pfinal)
  class(ret) <- c("classification_prediction", "prediction", "list")
  ret
  
}


#' Create a Model Fit Object
#'
#' Constructs a model fit object, representing the result of a single model fit to a chunk of data. The object contains information about the model, response variable, model fit, problem type, model parameters, voxel indices, and an optional feature mask.
#'
#' @param model The model specification object from MVPAModels.
#' @param y The response variable (predictand).
#' @param fit The fitted model.
#' @param model_type The problem type, either "classification" or "regression" (default). Must be one of the provided options.
#' @param param The model parameters.
#' @param vox_ind The voxel indices indicating the data coordinates.
#' @param feature_mask An optional logical mask indicating the selected subset of columns (features).
#'
#' @return An object of class \code{model_fit}, containing the model, response variable, fitted model, problem type, model parameters, voxel indices, and optional feature mask. The object is also assigned a class based on the problem type: \code{class_model_fit} for classification or \code{regression_model_fit} for regression.
#'
#' @keywords internal
#' @noRd
model_fit <- function(model, y, fit, model_type=c("classification", "regression"), param, vox_ind, feature_mask=NULL) {
  model_type=match.arg(model_type)
  
  ret <- list(
    model=model,
    y=y,
    fit=fit,
    model_type=model_type,
    param=param,
    vox_ind=vox_ind,
    feature_mask=feature_mask)
  
  if (model_type == "classification") {
    class(ret) <- c("class_model_fit", "model_fit")
  } else {
    class(ret) <- c("regression_model_fit", "model_fit")
  }
  ret
}

#' Create a Weighted Consensus Model
#'
#' Constructs a weighted consensus model formed as a weighted average of a set of models. The consensus model combines the input models according to their respective weights.
#'
#' @param fits A list of model fits to be combined.
#' @param names An optional list of names, one per model fit (default: numeric indices).
#' @param weights A vector of weights, one per model fit, that sum up to 1 (default: equal weights for all models).
#'
#' @return An object of class \code{weighted_model}, containing the list of model fits, their names, and the assigned weights. The object is also assigned a class `list`.
#'
#' @examples
#' # Create two sample model fits
#' fit1 <- list(model = "model1", y = c(0, 1), fit = "fit1")
#' fit2 <- list(model = "model2", y = c(1, 0), fit = "fit2")
#'
#' # Combine the model fits into a weighted consensus model
#' w_model <- weighted_model(fits = list(fit1, fit2), names = c("model1", "model2"), weights = c(0.6, 0.4))
#'
#' @keywords internal
#' @noRd
weighted_model <- function(fits, names=1:length(fits), weights=rep(1/length(fits), length(fits))) {
  stopifnot(length(weights) == length(fits))
  ret <- fits
  names(ret) <- names
  attr(ret, "weights") <- weights
  class(ret) <- c("weighted_model", "list")
  ret
}

#' a list of model fits
#' 
#' @param fits a list of fits
#' @param names the names of the fits
#' @noRd
#' @keywords internal
list_model <- function(fits, names=1:length(fits)) {
  stopifnot(is.list(fits))
  ret <- fits
  names(ret) <- names
  class(ret) <- c("list_model", "list")
  ret
}

#' @export
#' @method predict weighted_model
#' @importFrom stats predict
predict.weighted_model <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }

  preds <- lapply(object, function(fit) predict(fit, newdata, ...))
  merge_predictions(preds[[1]], preds[2:length(preds)], attr(object, "weights"))
  
}

#' @export
#' @method predict list_model
#' @importFrom stats predict
predict.list_model <- function(object, newdata=NULL,...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  res <- lapply(object, function(fit) {
    predict(fit, newdata,...)
  })
  
}



#' @rdname train_model
#' @param obj An object of class \code{mvpa_model}, specifying the MVPA problem.
#' @param train_dat Training data, an instance of class \code{ROIVolume} or \code{ROISurface}.
#' @param y The dependent variable (response variable), either a numeric vector or a factor.
#' @param indices The spatial indices associated with each column.
#' @param wts Optional class weights (if the underlying model supports it).
#' @param ... Additional arguments passed to other methods.
#' @return A model fit object containing the trained model, its fit, the model type (classification or regression), the best tuning parameters, the voxel indices, and the feature mask.
#' @method train_model mvpa_model
train_model.mvpa_model <- function(obj, train_dat, y, indices, wts=NULL, ...) {
  
  tryCatch({
    futile.logger::flog.debug("Starting train_model with data dimensions: %s", 
                             paste(dim(train_dat), collapse=" x "))
    futile.logger::flog.debug("Response variable levels: %s", 
                             paste(levels(y), collapse=", "))
    
    param <- tune_grid(obj, train_dat, y, len=1)
    futile.logger::flog.debug("Tuning grid parameters: %s", 
                             paste(names(param), collapse=", "))

    if (is.character(y)) {
      y <- as.factor(y)
    }
    
    ## columns that have zero variance
    nzero <- nonzeroVarianceColumns2(train_dat)
    futile.logger::flog.debug("Non-zero variance columns: %d", sum(nzero))
    
    ## columns with NAs
    nacols <- na_cols(train_dat)
    futile.logger::flog.debug("NA columns: %d", sum(nacols))
    
    ## duplicated columns
    dup <- !duplicated(t(train_dat))
    futile.logger::flog.debug("Non-duplicate columns: %d", sum(dup))
    
    ## invalid columns
    nzero <- nzero & dup & !nacols
    futile.logger::flog.debug("Valid columns after filtering: %d", sum(nzero))
    
    if (length(nzero) == 0 || sum(nzero,na.rm=TRUE) < 2) {
      stop(sprintf("training data must have more than one valid feature (found %d)", 
                  sum(nzero,na.rm=TRUE)))
    }
    
    ## feature selection and variable screening
    feature_mask <- if (!is.null(obj$feature_selector)) {
      nz <- which(nzero)
      fsel <- select_features(obj, train_dat[,nz], y)
      mask <- logical(ncol(train_dat))
      mask[nz[fsel]] <- TRUE
      mask
    } else {
      nzero
    }
    
    futile.logger::flog.debug("Features selected: %d", sum(feature_mask))
    
    if (sum(feature_mask) < 2) {
      stop("train_model: training data must have more than one valid feature after feature selection")
    }
    
    train_dat <- train_dat[,feature_mask,drop=FALSE]
    
    ## parameter_tuning
    best_param <- if (!is.vector(param) && !is.null(nrow(param)) && nrow(param) > 1) {
      bp <- tune_model(obj, train_dat, y, wts, param, obj$tune_reps)
      futile.logger::flog.debug("Best tuning parameters: %s", 
                               paste(capture.output(print(bp)), collapse="\n"))
      bp
    } else {
      param
    }
    
    mtype <- if (is.factor(y)) {
      "classification"
    } else if (is.numeric(y)) {
      "regression"
    } else {
      stop("'y' must be a numeric vector or factor")
    }
    
    futile.logger::flog.debug("Fitting model of type: %s", mtype)
    fit <- fit_model(obj, train_dat, y, wts=wts, param=best_param, classProbs=TRUE)
    model_fit(obj$model, y, fit, mtype, best_param, indices, feature_mask)
    
  }, error = function(e) {
    futile.logger::flog.error("train_model failed: %s", e$message)
    futile.logger::flog.debug("Data dimensions: %s", paste(dim(train_dat), collapse=" x "))
    futile.logger::flog.debug("Response levels: %s", paste(levels(y), collapse=", "))
    stop(e$message)  # Re-throw the error after logging
  })
}


