
#' binary_classification_result
#' 
#' create an \code{binary_classification_result} instance
#' 
#' @param observed the observed classes
#' @param predicted the predicted classes
#' @param probs the predicted probabilities
#' @export
binary_classification_result <- function(observed, predicted, probs, testDesign, predictor=NULL) {
  ret <- list(
    observed=observed,
    predicted=predicted,
    probs=as.matrix(probs),
    testDesign=testDesign,
    predictor=predictor
  )
  
  class(ret) <- c("binary_classification_result", "classification_result", "list")
  ret
}




sub_result.multiway_classification_result <- function(x, indices) {
  ret <- list(
    observed=x$observed[indices],
    predicted=x$predicted[indices],
    probs=as.matrix(x$probs)[indices,],
    testDesign=x$testDesign[indices,],
    predictor=x$predictor)
  
  class(ret) <- c("multiway_classification_result", "classification_result", "list")
  ret
}

sub_result.binary_classification_result <- function(x, indices) {
  ret <- list(
    observed=x$observed[indices],
    predicted=x$predicted[indices],
    probs=as.matrix(x$probs)[indices,],
    testDesign=x$testDesign[indices,],
    predictor=x$predictor)
  
  class(ret) <- c("binary_classification_result", "classification_result", "list")
  ret
}


#' create an \code{multiway_classification_result} instance
#' 
#' @param observed
#' @param predicted
#' @param probs
#' @export
multiway_classification_result <- function(observed, predicted, probs,  testDesign=NULL, predictor=NULL) {
  ret <- list(
    observed=observed,
    predicted=predicted,
    probs=as.matrix(probs),
    testDesign=testDesign,
    predictor=predictor)
  
  class(ret) <- c("multiway_classification_result", "classification_result", "list")
  ret
}

#' create an \code{regression_result} instance
#' 
#' @param observed the observed values
#' @param predicted the predicted values
#' @param testDesign
#' @param predictor
#' @export
regression_result <- function(observed, predicted, testDesign=NULL, predictor=NULL) {
  ret <- list(
    observed=observed,
    predicted=predicted,
    testDesign=testDesign,
    predictor=predictor)
  class(ret) <- c("regression_result", "classification_result", "list")
  ret
}



#' @export
classification_result <- function(observed, predicted, probs, testDesign=NULL,predictor=NULL) {
  if (is.numeric(observed)) {
    regression_result(observed, predicted, testDesign, predictor)
  } else if (length(levels(as.factor(observed))) == 2) {
    binary_classification_result(observed, predicted, probs,  testDesign, predictor)
  } else if (length(levels(as.factor(observed))) > 2) {
    multiway_classification_result(observed,predicted, probs, testDesign, predictor)
  } else {
    stop("observed data must be a factor with 2 or more levels")
  }
}
