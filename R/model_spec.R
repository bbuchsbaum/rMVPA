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
fit_model.model_spec <- function(obj, x, y, wts, param, lev, last, classProbs, ...) {
   obj$model$fit(x,y,wts,param,lev, last, classProbs, ...)
}

crossval_samples.model_spec <- function(obj) { crossval_samples(obj$crossval) }


#tune_model.model_spec <- function(obj, roi_train, Ytrain, blocking_var=NULL) {
#}
  

#' @param obj an instance of class \code{model_spec}
#' @param roi_train training data, and instance of class \code{ROIVolume} or \code{ROISurface}
#' @param Ytrain the dependent variable
#' @param wts
#' @export
train_model.model_spec <- function(obj, roi_train, Ytrain, param=NULL, wts=NULL) {
 
  if (!is.null(obj$feature_selector)) {
    feature_mask <- select_features(obj, roi_train, Ytrain)
    roi_train <- roi_train[feature_mask]
  } else {
    feature_mask <- rep(TRUE, length(roi_train))
  }

  if (is.null(param)) {
    param <- tune_grid(obj)[1,]
  }
  
  fit <- fit_model(obj, values(ROI), Ytrain, wts=wtsL, param=param, lev=levels(Ytrain), classProbs=TRUE)
  model_fit(obj$model, fit, param, indices(roi_train), feature_mask)
}

#' @export
predict.model_fit <- function(x, newdata, sub_indices) {
  mat <- if (is.matrix(newdata)) {
    newdata
  } else if (inherits(newdata, "BrainVector") || inherits(newdata, "BrainSurfaceVector")) {
    series(newdata, model_fit$vox_ind)
  } 
  
  if (!is.null(sub_indices)) {
    assert_that(is.vector(sub_indices))
    mat <- mat[sub_indices,,drop=FALSE]
  }
  
  probs <- x$model$prob(x$modelFit, mat[, x$featureMask,drop=FALSE]) 
  cpred <- max.col(probs)
  cpred <- colnames(probs)[cpred]
  list(class=cpred, probs=probs)
}


#' the result of a single model fit to a chunk of data
#' 
#' @param model the caret-style model object
#' @param fit the model fit 
#' @param param the model parameters
#' @param vox_ind the the voxel indices indicating the data coordinates
#' @param feature_mask a logical mask indicating the selected subset of columns
#' @export
model_fit <- function(model, fit, param, vox_ind, feature_mask=NULL) {
  ret <- list(
    model=model,
    fit=fit,
    param=param,
    vox_ind=vox_ind,
    feature_mask=feature_mask)
  
  class(ret) <- "model_fit"
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
