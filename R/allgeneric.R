

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
#' @export
performance <- function(x,...) {
  UseMethod("performance")
}


#' compute_performance
#' 
#' Delegate calculation of performance metrics
#' 
compute_performance <- function(x, result) {
  UseMethod("compute_performance")
}

#' merge_results
#' 
#' merge two classification results
#' 
#' @param x the first result
#' @param y the second result
#' @export
merge_results <- function(x, y, ...) {
  UseMethod("merge_results")
}


#' get_samples
#' 
#' @export
get_samples <- function(obj, voxiter) {
  UseMethod("get_samples")
}

#' data_sample
#' 
#' @export
data_sample <- function(obj, vox) {
  UseMethod("data_sample")
}


#' extract_sample
#' 
#' @export
extract_sample <- function(obj,...) {
  UseMethod("extract_sample")
}


#' as_roi
#' 
#' @export
as_roi <- function(obj,...) {
  UseMethod("as_roi")
}


#' get_searchlight
#' 
get_searchlight <- function(obj, ...) {
  UseMethod("get_searchlight")
}


#' wrap_output
#' 
wrap_output <- function(obj, vals, ...) {
  UseMethod("wrap_output")
}



