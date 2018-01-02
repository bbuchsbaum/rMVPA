#' create a \code{classification_result} instance
#' 
#' @param observed the observed values.
#' @param predicted the predicted values.
#' @param probs a \code{matrix} of predicted probabilities, one column per level.
#' @param test_design an optional design for the test data.
#' @param predictor an optional predictor object.
#' @rdname classification_result
#' @export
classification_result <- function(observed, predicted, probs, test_design=NULL,predictor=NULL) {
  if (is.numeric(observed)) {
    regression_result(observed, predicted, test_design, predictor)
  } else if (length(levels(as.factor(observed))) == 2) {
    binary_classification_result(observed, predicted, probs,  test_design, predictor)
  } else if (length(levels(as.factor(observed))) > 2) {
    multiway_classification_result(observed,predicted, probs, test_design, predictor)
  } else {
    stop("observed data must be a factor with 2 or more levels")
  }
}


#' @inheritParams classification_result
#' @rdname classification_result
#' @export
binary_classification_result <- function(observed, predicted, probs, test_design=NULL, predictor=NULL) {
  ret <- list(
    observed=observed,
    predicted=predicted,
    probs=as.matrix(probs),
    test_design=test_design,
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
    test_design=x$test_design[indices,],
    predictor=x$predictor)
  
  class(ret) <- c("multiway_classification_result", "classification_result", "list")
  ret
}

sub_result.binary_classification_result <- function(x, indices) {
  ret <- list(
    observed=x$observed[indices],
    predicted=x$predicted[indices],
    probs=as.matrix(x$probs)[indices,],
    test_design=x$test_design[indices,],
    predictor=x$predictor)
  
  class(ret) <- c("binary_classification_result", "classification_result", "list")
  ret
}


 
#' @inheritParams classification_result
#' @rdname classification_result
#' @export
multiway_classification_result <- function(observed, predicted, probs,  test_design=NULL, predictor=NULL) {
  ret <- list(
    observed=observed,
    predicted=predicted,
    probs=as.matrix(probs),
    test_design=test_design,
    predictor=predictor)
  
  class(ret) <- c("multiway_classification_result", "classification_result", "list")
  ret
}

 
#' @inheritParams classification_result
#' @rdname classification_result
#' @export
regression_result <- function(observed, predicted, test_design=NULL, predictor=NULL) {
  ret <- list(
    observed=observed,
    predicted=predicted,
    test_design=test_design,
    predictor=predictor)
  class(ret) <- c("regression_result", "classification_result", "list")
  ret
}

