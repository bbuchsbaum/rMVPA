
#' @keywords internal
requireNamespaceQuietStop <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop(paste('package',package,'is required'), call. = FALSE)
}

#' @keywords internal
mclass_summary <- function (data, lev = NULL, model = NULL) {
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  has_class_probs <- all(lev %in% colnames(data))
  
  if (has_class_probs) {
    requireNamespaceQuietStop("ModelMetrics")
    prob_stats <- lapply(levels(data[, "pred"]), function(x) {
      obs <- ifelse(data[, "obs"] == x, 1, 0)
      prob <- data[, x]
      AUCs <- try(ModelMetrics::auc(obs, data[, x]), silent = TRUE)
      return(AUCs)
    })
    roc <- mean(unlist(prob_stats))
  } else {
    stop("Cannot compute AUC. Class probabilities unavailable for model: ", model)
  }
  
  c(AUC=roc)
}



#' @keywords internal
load_caret_libs <- function(x) {
  for (lib in x$model$library) {
    library(lib, character.only = TRUE)
  }
}


#' @keywords internal
get_control <- function(y, nreps) {
  if (is.factor(y) && length(levels(y)) == 2) {
    ctrl <- caret::trainControl("boot", number=nreps, verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE, trim=TRUE, summaryFunction=caret::twoClassSummary)
    metric <- "ROC"
  } else if (is.factor(y) && length(levels(y)) > 2) {
    ctrl <- caret::trainControl("boot", number=nreps, verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE, trim=TRUE, summaryFunction=mclass_summary)
    metric <- "AUC"
  } else {
    ctrl <- caret::trainControl("boot", number=nreps, verboseIter=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE, trim=TRUE)
    metric = "RMSE"
  }
  
  list(ctrl=ctrl, metric=metric)
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
tune_model <- function(mspec, x, y, wts, param, nreps=10) {
  ctrl <- get_control(y, nreps)
  cfit <-caret::train(as.data.frame(x), y, method=mspec$model, weights=wts, metric=ctrl$metric, trControl=ctrl$ctrl, tuneGrid=param)
  cfit$bestTune
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
fit_model.mvpa_model <- function(obj, x, y, wts, param, classProbs, ...) {
  fit <- obj$model$fit(x,y,wts=wts,param=param,lev=levels(y), classProbs=classProbs, ...)
  ##fit$obsLevels <- levels(y)
  attr(fit, "obsLevels") <- levels(y)
  if (is.factor(y)) {
    #fit$problemType <- "Classification"
    attr(fit, "problemType") <- "Classification"
  } else {
    #fit$problemType <- "Regression"
    attr(fit, "problemType") <- "Regression"
  }
  fit
}

#'
#' This function predicts class labels and probabilities for new data using a fitted model.
#'
#' @param object A fitted model object of class \code{class_model_fit}.
#' @param newdata New data to predict on, either as a \code{matrix} or a \code{NeuroVec} or \code{NeuroSurfaceVector} object.
#' @param sub_indices The subset of row indices to compute predictions on (optional).
#' @param ... Additional arguments to be passed to the underlying prediction function.
#' @return A list containing class predictions and probabilities with class attributes "classification_prediction", "prediction", and "list".
predict.class_model_fit <- function(object, newdata, sub_indices=NULL,...) {
  
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

  probs <- object$model$prob(object$fit,mat) 
  colnames(probs) <- levels(object$y)
  cpred <- max.col(probs)
  cpred <- levels(object$y)[cpred]
  ret <- list(class=cpred, probs=probs)
  class(ret) <- c("classification_prediction", "prediction", "list")
  ret
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
#' @examples
#' # Assuming `fitted_model` is a fitted model object of class `regression_model_fit`
#' new_data <- iris_dataset$test_data
#' predictions <- predict(fitted_model, new_data)
#' @export
predict.regression_model_fit <- function(object, newdata, sub_indices=NULL,...) {
  
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
  
  ret <- list(preds=object$model$predict(object$fit,mat))
  class(ret) <- c("regression_prediction", "prediction", "list")
  ret
}


#' @export
merge_predictions.regression_prediction <- function(obj1, rest, weights=rep(1,length(rest)+1)/(length(rest)+1)) {
  allobj <- c(obj1, rest)
  assert_that(all(sapply(allobj, function(obj) inherits(obj, "regression_prediction"))))
  
  #preds <- lapply(1:length(allobj), function(i) {
  #  predict(allobj[[i]], newdata, ...)$pred * weights[i]
  #})
  
  preds <- lapply(1:length(allobj), function(i) {
    allobj[[i]]$pred * weights[i]
  })
  
  final_pred <- rowMeans(do.call(cbind, preds))
  ret <- list(preds=final_pred)
  class(ret) <- c("regression_prediction", "prediction", "list")
  ret
}


#' @export
#' @method merge_predictions classification_prediction
merge_predictions.classification_prediction <- function(obj1, rest, weights=rep(1,length(rest)+1)/(length(rest)+1)) {
  allobj <- vector(mode="list", length(rest)+1)
  allobj[[1]] <- obj1
  allobj[2:length(allobj)] <- rest
  
  #allobj <- c(obj1, rest)
  assert_that(all(sapply(allobj, function(obj) inherits(obj, "classification_prediction"))))
  
  #preds <- lapply(1:length(allobj), function(i) {
  #  predict(allobj[[i]], newdata, ...)$prob * weights[i]
  #})
  
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
#' @param model The caret-style model object.
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
#' @export
list_model <- function(fits, names=1:length(fits)) {
  stopifnot(is.list(fits))
  ret <- fits
  names(ret) <- names
  class(ret) <- c("list_model", "list")
  ret
}

#' @export
#' @method predict weighted_model
predict.weighted_model <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }

  preds <- lapply(object, function(fit) predict(fit, newdata, ...))
  merge_predictions(preds[[1]], preds[2:length(preds)], attr(object, "weights"))
  
}

#' @export
#' @method predict list_model
predict.list_model <- function(object, newdata=NULL,...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  res <- lapply(object, function(fit) {
    predict(fit, newdata,...)
  })
  
}



#' Train an MVPA Model
#'
#' This function trains a Multi-Variate Pattern Analysis (MVPA) model on the provided data, taking care of feature selection, parameter tuning, and model fitting.
#'
#' @param obj An object of class \code{mvpa_model}, specifying the MVPA problem.
#' @param train_dat Training data, an instance of class \code{ROIVolume} or \code{ROISurface}.
#' @param y The dependent variable (response variable), either a numeric vector or a factor.
#' @param indices The spatial indices associated with each column.
#' @param param Optional tuning parameters (default: NULL). If not provided, the function will perform tuning using \code{tune_grid}.
#' @param wts Optional case weights (default: NULL).
#' @param tune_reps The number of bootstrap replications for parameter tuning (default: 10), only used when \code{param} is NULL.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A model fit object containing the trained model, its fit, the model type (classification or regression), the best tuning parameters, the voxel indices, and the feature mask.
train_model.mvpa_model <- function(obj, train_dat, y, indices, param=NULL, wts=NULL, tune_reps=10,...) {
  
  if (is.null(param)) {
    param <- tune_grid(obj, train_dat, y, len=1)
  }
  
  if (is.character(y)) {
    y <- as.factor(y)
  }
  
 
  ## columns that have zero variance
  nzero <- nonzeroVarianceColumns2(train_dat)
  ## columns with NAs
  nacols <- na_cols(train_dat)
  ## duplicated columns
  dup <- !duplicated(t(train_dat))
  ## invalid columns
  nzero <- nzero & dup & !nacols
  
  if (length(nzero) == 0 || sum(nzero,na.rm=TRUE) < 2) {
      stop("training data must have more than one valid feature")
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
  
  if (sum(feature_mask) < 2) {
    stop("train_model: training data must have more than one valid feature")
  }
  
  #print(feature_mask)
  
  train_dat <- train_dat[,feature_mask]
  
  ## parameter_tuning
  best_param <- if (!is.vector(param) && !is.null(nrow(param)) && nrow(param) > 1) {
    bp <- tune_model(obj, train_dat, y, wts, param, tune_reps)
    flog.info("best tuning parameter: ", bp, capture=TRUE)
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
  
  fit <- fit_model(obj, train_dat, y, wts=wts, param=best_param, classProbs=TRUE)
  model_fit(obj$model, y, fit, mtype, best_param, indices, feature_mask)
}


