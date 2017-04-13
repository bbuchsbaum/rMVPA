

#' select_features
#' 
#' performance feature selection on a training ROI.
#' 
#' @export
select_features <- function(obj, roi, Y, ...) {
  UseMethod("select_features")
}


#' train_model
#' 
#' train a classification or regression model
#' 
#' @export
train_model <- function(obj,...) {
  UseMethod("train_model")
}



#' y_train
#' 
#' extract the training labels/response
#' 
#' @export
y_train <- function(obj) {
  UseMethod("y_train")
}


#' y_test
#' 
#' extract the test labels/response.
#' 
#' @export
y_test <- function(obj) {
  UseMethod("y_test")
}

#' fit_model
#' 
#' fit a classification or regression model
#' 
fit_model <- function(obj, roi_x, y, wts, param, lev, last, classProbs, ...) {
  UseMethod("fit_model")
}


#' tune_grid
#' 
#' extract the parameter grid to optimize.
#' 
#' @export
tune_grid <- function(obj, x,y,len) {
  UseMethod("tune_grid")
}

#' has_test_set
#' 
#' @export
has_test_set <- function(obj) {
  UseMethod("has_test_set")
}


#' performance
#' 
#' Compute appropriate performance metrics such as accuracy/AUC/RMSE from a classification/regression result
#' 
#' @param x the result to evaluate performance of
#' @param ... extra args
#' @export
performance <- function(x,...) {
  UseMethod("performance")
}


#' compute_performance
#' 
#' Delegate calculation of performance metrics
#' 
#' @param obj the object
#' @param result the 'result' to evaluate
compute_performance <- function(obj, result) {
  UseMethod("compute_performance")
}

#' merge_results
#' 
#' merge two classification/regression results
#' 
#' @param x the first result
#' @param y the second result
#' @param extra args
#' @export
merge_results <- function(x, y, ...) {
  UseMethod("merge_results")
}


#' get_samples
#' 
#' @param obj the object
#' @param voxiter the list of voxel/index sets
#' @export
get_samples <- function(obj, voxiter) {
  UseMethod("get_samples")
}

#' data_sample
#' 
#' @param obj the object
#' @param vox the voxel indices
#' @export
data_sample <- function(obj, vox) {
  UseMethod("data_sample")
}


#' extract_sample
#' 
#' @param obj the object
#' @param ... extra args
extract_sample <- function(obj,...) {
  UseMethod("extract_sample")
}


#' as_roi
#' 
#' convert to an \code{ROIVolume} or \code{ROISurface} object
#' 
#' @param obj the object to convert
#' @param ... extra args
#' @export
as_roi <- function(obj,...) {
  UseMethod("as_roi")
}


#' get_searchlight
#' 
#' get a searchlight object appropriate for a given input dataset (e.g. volumetric or surface).
#' 
#' @param obj the object
#' @param ... extra args
#' @export
get_searchlight <- function(obj, ...) {
  UseMethod("get_searchlight")
}


#' wrap_output
#' 
#' @param obj the object
#' @param vals the values to wrap
#' @param ... extra args
wrap_output <- function(obj, vals, ...) {
  UseMethod("wrap_output")
}

#' merge_predictions
#' 
#' combine predictions from several models applied to the same test set.
#' 
#' @param obj1 the first object
#' @param rest the rest of the objects
#' @param ... extra args
merge_predictions <- function(obj1, rest, ...) {
  UseMethod("merge_predictions")
}


#' sub_result
#' 
#' exctract a row-wise subset of a classification/regression result object.
#' 
#' @param x the result object to subset
#' @param indices the row indices
sub_result <- function(x, indices) {
  UseMethod("sub_result")
}



