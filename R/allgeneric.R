


#' @export
select_features <- function(obj, roi, Y, ...) {
  UseMethod("select_features")
}

#' @export
train_model <- function(obj,...) {
  UseMethod("train_model")
}

#' @export
y_train <- function(obj) {
  UseMethod("y_train")
}

#' @export
y_test <- function(obj) {
  UseMethod("y_test")
}


#' @export
fit_model <- function(obj, roi_x, y, wts, param, lev, last, classProbs, ...) {
  UseMethod("fit_model")
}

#' @export
tune_grid <- function(obj, x,y,len) {
  UseMethod("tune_grid")
}


#' @export
has_test_set <- function(obj) {
  UseMethod("has_test_set")
}


#' performance
#' 
#' Compute performance metrics from a classiifcation result
#' 
#' @param x the result to evaluate performance of
#' @export
performance <- function(x,...) {
  UseMethod("performance")
}

#' merge_results
#' 
#' merge two classiifcation results
#' 
#' @param x the first result
#' @param y the second result
#' @export
merge_results <- function(x, y, ...) {
  UseMethod("merge_results")
}

#' @export
get_samples <- function(obj, voxiter) {
  UseMethod("get_samples")
}


#' @export
data_sample <- function(obj, vox) {
  UseMethod("data_sample")
}

#' @export
extract_sample <- function(obj,...) {
  UseMethod("extract_sample")
}

#' @export
as_roi <- function(obj,...) {
  UseMethod("as_roi")
}



