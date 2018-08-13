
mclass_summary <- function (data, lev = NULL, model = NULL) {
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  has_class_probs <- all(lev %in% colnames(data))
  
  if (has_class_probs) {
    caret::requireNamespaceQuietStop("ModelMetrics")
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




load_caret_libs <- function(x) {
  for (lib in x$model$library) {
    library(lib, character.only = TRUE)
  }
}



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


#' tune_model
#' 
#' find the best tuning parameters for a model specification
#' 
#' @param mspec the model specification which derives from class \code{mvpa_model}
#' @param x the training data matrix
#' @param y the response vector
#' @param wts optional class weights (if underlying model supports it)
#' @param param the tuning grid, should be a \code{data.frame} where parameter names are indicated by column names.
#' @param nreps the number of bootstrap replications
tune_model <- function(mspec, x, y, wts, param, nreps=10) {
  ctrl <- get_control(y, nreps)
  cfit <-caret::train(as.data.frame(x), y, method=mspec$model, weights=wts, metric=ctrl$metric, trControl=ctrl$ctrl, tuneGrid=param)
  cfit$bestTune
}


fit_model.mvpa_model <- function(obj, x, y, wts, param, classProbs, ...) {
  obj$model$fit(x,y,wts=wts,param=param,classProbs=classProbs, ...)
}



#' @export
#' @param sub_indices the subset of row indices to compute predictions on
predict.class_model_fit <- function(object, newdata, sub_indices=NULL,...) {
  
  mat <- if (inherits(newdata, "BrainVector") || inherits(newdata, "BrainSurfaceVector")) {
    series(newdata, x$fit$vox_ind)
  } else {
    newdata
  }
  
  if (!is.null(sub_indices)) {
    assert_that(is.vector(sub_indices))
    mat <- mat[sub_indices,,drop=FALSE]
  }
  
  if (!is.null(x$feature_mask)) {
    mat <- mat[, x$feature_mask,drop=FALSE]
  }
  
  probs <- x$model$prob(x$fit,mat) 
  colnames(probs) <- levels(x$y)
  cpred <- max.col(probs)
  cpred <- levels(x$y)[cpred]
  ret <- list(class=cpred, probs=probs)
  class(ret) <- c("classification_prediction", "prediction", "list")
  ret
}

#' @export
predict.regression_model_fit <- function(object, newdata, sub_indices=NULL) {
  
  mat <- if (inherits(newdata, "BrainVector") || inherits(newdata, "BrainSurfaceVector")) {
    series(newdata, x$fit$vox_ind)
  } else {
    newdata
  }
  
  if (!is.null(sub_indices)) {
    assert_that(is.vector(sub_indices))
    mat <- mat[sub_indices,,drop=FALSE]
  }
  
  
  if (!is.null(x$feature_mask)) {
    mat <- mat[, x$feature_mask,drop=FALSE]
  }
  
  ret <- list(preds=x$model$predict(x$fit,mat))
  class(ret) <- c("regression_prediction", "prediction", "list")
  ret
}


#' @export
merge_predictions.regression_prediction <- function(obj1, rest, weights=rep(1,length(rest)+1)/(length(rest)+1)) {
  allobj <- c(obj1, rest)
  assert_that(all(sapply(allobj, function(obj) inherits(obj, "regression_prediction"))))
  
  preds <- lapply(1:length(allobj), function(i) {
    predict(allobj[[i]], newdata, ...)$pred * weights[i]
  })
  
  final_pred <- rowMeans(do.call(cbind, preds))
  ret <- list(preds=final_pred)
  class(ret) <- c("regression_prediction", "prediction", "list")
  ret
}


#' @export
merge_predictions.classification_prediction <- function(obj1, rest, weights=rep(1,length(rest)+1)/(length(rest)+1)) {
  allobj <- c(obj1, rest)
  assert_that(all(sapply(allobj, function(obj) inherits(obj, "classification_prediction"))))
  
  
  preds <- lapply(1:length(allobj), function(i) {
    predict(allobj[[i]], newdata, ...)$prob * weights[i]
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


#' the result of a single model fit to a chunk of data
#' 
#' @param model the caret-style model object
#' @param y the predictand
#' @param fit the model fit 
#' @param model_type the problem type: classification or regression
#' @param param the model parameters
#' @param vox_ind the the voxel indices indicating the data coordinates
#' @param feature_mask a logical mask indicating the selected subset of columns
#' @export
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

#' weighted_model
#' 
#' a consensus model formed as a weighted average of a set of models
#' 
#' @param fits a list of model fits
#' @param names a list of names, one per model fit
#' @param weights a vector of weights, one per model fit
#' @export
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
predict.weighted_model <- function(object, newdata=NULL, ...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }

  preds <- lapply(x$fits, function(fit) predict(fit, newdata, ...))
  merge_predictions(preds[[1]], preds[2:length(preds)], attr(x, "weights"))
  
}

#' @export
predict.list_model <- function(object, newdata=NULL,...) {
  if (is.null(newdata)) {
    stop("newdata cannot be null")
  }
  
  res <- lapply(x, function(fit) {
    predict(fit, newdata,...)
  })
  
}






#' train_model
#' 
#' @param train_dat training data, and instance of class \code{ROIVolume} or \code{ROISurface}
#' @param y the dependent variable
#' @param indices the spatial indices associated with each column
#' @param param optional tuning parameters
#' @param wts optional case weights
#' @param tune_reps the number of bootstrap replications for parameter tuning (only used when param is not \code{NULL})
#' @export
#' @describeIn train_model train an mvpa_model
train_model.mvpa_model <- function(obj, train_dat, y, indices, param=NULL, wts=NULL, tune_reps=10,...) {
  
  if (is.null(param)) {
    param <- tune_grid(obj, train_dat, y, len=1)
  }
  
  if (is.character(y)) {
    y <- as.factor(y)
  }
  
  nzero <- nonzeroVarianceColumns2(train_dat)
  dup <- !duplicated(t(train_dat))
  
  nzero <- nzero & dup
  
  if (sum(nzero) < 2) {
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
    stop("training data must have more than one valid feature")
  }
  
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


