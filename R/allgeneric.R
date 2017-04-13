

#' select_features
#' 
#' @export
select_features <- function(obj, roi, Y, ...) {
  UseMethod("select_features")
}


#' train_model
#' 
#' @export
train_model <- function(obj,...) {
  UseMethod("train_model")
}



#' y_train
#' 
#' @export
y_train <- function(obj) {
  UseMethod("y_train")
}


#' y_test
#' 
#' @export
y_test <- function(obj) {
  UseMethod("y_test")
}

#' fit_model
#' 
#' @export
fit_model <- function(obj, roi_x, y, wts, param, lev, last, classProbs, ...) {
  UseMethod("fit_model")
}


#' tune_grid
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
#' Calculate performance metrics from a classification/regression result
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
#' merge two classification results
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
#' @param obj the object to convert
#' @param ... extra args
#' @export
as_roi <- function(obj,...) {
  UseMethod("as_roi")
}


#' get_searchlight
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
#' @param obj1 the first object
#' @param rest the rest of the objects
#' @param ... extra args
merge_predictions <- function(obj1, rest, ...) {
  UseMethod("merge_predictions")
}


#' sub_result
#' @param x the result object to subset
#' @param indices the row indices
sub_result <- function(x, indices) {
  UseMethod("sub_result")
}



