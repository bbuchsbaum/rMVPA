
mclass_summary <- function (data, lev = NULL, model = NULL) {
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"]))) 
    stop("levels of observed and predicted data do not match")
  has_class_probs <- all(lev %in% colnames(data))
  
  if (has_class_probs) {
    caret:::requireNamespaceQuietStop("ModelMetrics")
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



#' @export
load_libs.caret_model_wrapper <- function(x) {
  for (lib in x$model$library) {
    library(lib, character.only = TRUE)
  }
}


#' @export
tune_grid.model_spec <- function(obj, x,y,len) {
  if (is.null(obj$tune_grid)) {
    obj$model$grid(x,y,len)
  } else {
    obj$tune_grid
  }
}


#' @export
select_features.model_spec <- function(obj, roi, Y) {
  if (!is.null(obj$feature_selector)) {
    selectFeatures(obj$feature_selector, roi, Y)
  } else {
    rep(TRUE, length(roi))
  }
}


#' @export
crossval_samples.model_spec <- function(obj) { crossval_samples(obj$crossval) }

get_control <- function(y, nreps) {
  if (is.factor(y) && length(levels(y)) == 2) {
    ctrl <- caret::trainControl("boot", number=nreps, verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE, trim=TRUE, summaryFunction=caret::twoClassSummary)
    metric <- "AUC"
  } else if (is.factor(y) && length(levels(y)) > 2) {
    ctrl <- caret::trainControl("boot", number=nreps, verboseIter=TRUE, classProbs=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE, trim=TRUE, summaryFunction=mclass_summary)
    metric <- "AUC"
  } else {
    ctrl <- caret::trainControl("boot", number=nreps, verboseIter=TRUE, returnData=FALSE, returnResamp="none",allowParallel=FALSE, trim=TRUE)
    metric = "RMSE"
  }
  
  list(ctrl=ctrl, metric=metric)
}

tune_model <- function(mspec, x, y, wts, param, nreps=2) {
  ctrl <- get_control(y, nreps)
  cfit <-caret::train(as.data.frame(x), y, method=mspec$model, weights=wts, metric=ctrl$metric, trControl=ctrl$ctrl, tuneGrid=param)
  cfit$bestTune
}

#' @export
fit_model.model_spec <- function(obj, x, y, wts, param, classProbs, ...) {
  obj$model$fit(x,y,wts=wts,param=param,classProbs=classProbs, ...)
}



#' train_model
#' 
#' @param obj an instance of class \code{model_spec}
#' @param train_dat training data, and instance of class \code{ROIVolume} or \code{ROISurface}
#' @param y the dependent variable
#' @param indices the spatial indices associated with each column
#' @param param
#' @param wts
#' @export
train_model.model_spec <- function(obj, train_dat, y, indices, param=NULL, wts=NULL) {
  
  if (is.null(param)) {
    param <- tune_grid(obj)
  }
  
  if (is.character(y)) {
    y <- as.factor(y)
  }
  
  
  ## feature selection
  if (!is.null(obj$feature_selector)) {
    feature_mask <- select_features(obj, train_dat, y)
    train_dat <- train_dat[,feature_mask]
  } else {
    feature_mask <- rep(TRUE, ncol(train_dat))
  }
  
  
  ## parameter_tuning
  best_param <- if (!is.vector(param) && !is.null(nrow(param)) && nrow(param) > 1) {
    tune_model(obj, train_dat[, feature_mask], y, wts, param)
  } else {
    param
  }
  
  print(best_param)
  
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

#' @export
predict.class_model_fit <- function(x, newdata, sub_indices=NULL) {
  
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
  names(probs) <- levels(x$y)
  cpred <- max.col(probs)
  cpred <- levels(x$y)[cpred]
  list(class=cpred, probs=probs)
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

#' @param model
#' @param model_type
#' @param crossval
#' @param feature_selector
#' @param tune_grid
#' @param custom_performance
#' @export
model_spec <- function(model, model_type=c("classification", "regression"), 
                       crossval, 
                       feature_selector=NULL, tune_grid=NULL, 
                       custom_performance=NULL) {
  
  if (!is.null(custom_performance)) {
    assert_that(is.function(custom_performance)) 
  }
  
  model_type <- match.arg(model_type)
  
  ret <- list(model=model,
              model_type=model_type,
              model_name=model$label,
              tune_grid=tune_grid,
              feature_selector=feature_selector,
              crossval=crossval,
              custom_performance=custom_performance)
  
  class(ret) <- "model_spec"
  ret
  
}

