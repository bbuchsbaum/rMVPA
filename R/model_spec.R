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
fit_model.model_spec <- function(obj, x, y, wts, param, classProbs, ...) {
   obj$model$fit(x,y,wts=wts,param=param,classProbs=classProbs, ...)
}

crossval_samples.model_spec <- function(obj) { crossval_samples(obj$crossval) }


#tune_model.model_spec <- function(obj, roi_train, Ytrain, blocking_var=NULL) {
#}
  

#' @param obj an instance of class \code{model_spec}
#' @param train_dat training data, and instance of class \code{ROIVolume} or \code{ROISurface}
#' @param y the dependent variable
#' @param indices the spatial indices associated with each column
#' @param
#' @param wts
#' @export
train_model.model_spec <- function(obj, train_dat, y, indices, param=NULL, wts=NULL) {
 
  if (is.null(param)) {
    param <- tune_grid(obj)[1,]
  }
  
  if (is.character(y)) {
    y <- as.factor(y)
  }
  
  mtype <- if (is.factor(y)) {
    "classification"
  } else if (is.numeric(y)) {
    "regression"
  } else {
    stop("'y' must be a numeric vector or factor")
  }
  
  if (!is.null(obj$feature_selector)) {
    feature_mask <- select_features(obj, train_dat, y)
    train_dat <- train_dat[,feature_mask]
  } else {
    feature_mask <- rep(TRUE, ncol(train_dat))
  }
  
  
  fit <- fit_model(obj, train_dat, y, wts=wts, param=param, classProbs=TRUE)
  model_fit(obj$model, y, fit, mtype, param, indices, feature_mask)
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
#' @param crossval
#' @param feature_selector
#' @param tune_grid
#' @param custom_performance
#' @export
model_spec <- function(model, crossval, feature_selector=NULL, tune_grid=NULL, custom_performance=NULL) {
  
  if (!is.null(custom_performance)) {
    assert_that(is.function(custom_performance)) 
  }
  
  ret <- list(model=model,
       model_name=model$label,
       tune_grid=tune_grid,
       feature_selector=feature_selector,
       crossval=crossval,
       custom_performance=custom_performance)
  
  class(ret) <- "model_spec"
  ret
  
}


#mvpa_crossval <- function()
