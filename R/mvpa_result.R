#' Create a \code{classification_result} instance
#'
#' Constructs a classification result object based on the observed and predicted values,
#' as well as other optional parameters.
#'
#' @param observed A vector of observed or true values.
#' @param predicted A vector of predicted values.
#' @param probs A \code{matrix} of predicted probabilities, with one column per level.
#' @param testind The row indices of the test observations (optional).
#' @param test_design An optional design for the test data.
#' @param predictor An optional predictor object.
#'
#' @return A classification result object, which can be one of: \code{regression_result},
#'   \code{binary_classification_result}, or \code{multiway_classification_result}.
#'
#' @examples
#' # A vector of observed values
#' yobs <- factor(rep(letters[1:4], 5))
#'
#' # Predicted probabilities
#' probs <- data.frame(a = runif(1:20), b = runif(1:20), c = runif(1:20), d = runif(1:20))
#' probs <- sweep(probs, 1, rowSums(probs), "/")
#'
#' # Get the max probability per row and use this to determine the predicted class
#' maxcol <- max.col(probs)
#' predicted <- levels(yobs)[maxcol]
#'
#' # Construct a classification result
#' cres <- classification_result(yobs, predicted, probs)
#'
#' # Compute default performance measures (Accuracy, AUC)
#' performance(cres)
#' @export
classification_result <- function(observed, predicted, probs, testind=NULL, test_design=NULL,predictor=NULL) {
  
  
  if (is.numeric(observed)) {
    regression_result(observed, predicted, testind, test_design, predictor)
  } else if (length(levels(as.factor(observed))) == 2) {
    binary_classification_result(as.factor(observed), predicted, probs,  testind, test_design, predictor)
  } else if (length(levels(as.factor(observed))) > 2) {
    multiway_classification_result(as.factor(observed),predicted, probs, testind, test_design, predictor)
  } else {
    stop("observed data must be a factor with 2 or more levels")
  }
}

#' Classification results for binary outcome
#'
#' Constructs a binary classification result object based on the observed and predicted values,
#' as well as other optional parameters.
#'
#' @param observed A vector of observed or true values.
#' @param predicted A vector of predicted values.
#' @param probs A \code{matrix} of predicted probabilities, with one column per level.
#' @param testind The row indices of the test observations (optional).
#' @param test_design An optional design for the test data.
#' @param predictor An optional predictor object.
#'
#' @return A binary classification result object, with the class attribute set to "binary_classification_result".
#' @rdname classification_result
#' @export
binary_classification_result <- function(observed, predicted, probs, testind=NULL, test_design=NULL, predictor=NULL) {
  assertthat::assert_that(length(observed) == length(predicted))
  ret <- list(
    observed=observed,
    predicted=predicted,
    probs=as.matrix(probs),
    testind=testind,
    test_design=test_design,
    predictor=predictor
  )
  
  class(ret) <- c("binary_classification_result", "classification_result", "list")
  ret
}



#' @export
sub_result.multiway_classification_result <- function(x, indices) {
  ret <- list(
    observed=x$observed[indices],
    predicted=x$predicted[indices],
    probs=as.matrix(x$probs)[indices,],
    testind=x$testind[indices],
    test_design=x$test_design[indices,],
    predictor=x$predictor)
  
  class(ret) <- c("multiway_classification_result", "classification_result", "list")
  ret
}

#' @export
#' @param indices the set of indices used to subset results
#' @rdname sub_result
sub_result.binary_classification_result <- function(x, indices) {
  ret <- list(
    observed=x$observed[indices],
    predicted=x$predicted[indices],
    probs=as.matrix(x$probs)[indices,],
    testind=x$testind[indices],
    test_design=x$test_design[indices,],
    predictor=x$predictor)
  
  class(ret) <- c("binary_classification_result", "classification_result", "list")
  ret
}


 
#' @inheritParams classification_result
#' @rdname classification_result
#' @export
multiway_classification_result <- function(observed, predicted, probs,testind=NULL, test_design=NULL, predictor=NULL) {
  assertthat::assert_that(length(observed) == length(predicted))
  ret <- list(
    observed=observed,
    predicted=predicted,
    probs=as.matrix(probs),
    testind=testind,
    test_design=test_design,
    predictor=predictor)
  
  class(ret) <- c("multiway_classification_result", "classification_result", "list")
  ret
}

 
#' @inheritParams classification_result
#' @rdname classification_result
#' @export
regression_result <- function(observed, predicted, testind=NULL, test_design=NULL, predictor=NULL) {
  ret <- list(
    observed=observed,
    predicted=predicted,
    test_design=test_design,
    testind=testind,
    predictor=predictor)
  class(ret) <- c("regression_result", "classification_result", "list")
  ret
}

