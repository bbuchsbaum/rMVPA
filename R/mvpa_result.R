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
#' @family classification_result
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
#' @family classification_result
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



#' Subset Multiway Classification Result
#'
#' This function subsets a multiway classification result based on the provided indices.
#'
#' @param x An object of class \code{multiway_classification_result} containing the multiway classification results.
#' @param indices The set of indices used to subset the results.
#'
#' @return A \code{multiway_classification_result} object containing the subset of results specified by the indices.
#'
#' @export
#' @family sub_result
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

#' Subset Binary Classification Result
#'
#' This function subsets a binary classification result based on the provided indices.
#'
#' @param x An object of class \code{binary_classification_result} containing the binary classification results.
#' @param indices The set of indices used to subset the results.
#'
#' @return A \code{binary_classification_result} object containing the subset of results specified by the indices.
#'
#' @export
#' @family sub_result
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


 
#' Create a Multiway Classification Result Object
#'
#' This function creates a multiway classification result object containing the observed and predicted values, class probabilities, test design, test indices, and predictor.
#'
#' @param observed A vector of observed values.
#' @param predicted A vector of predicted values.
#' @param probs A matrix of class probabilities.
#' @param testind A vector of indices for the test data (optional).
#' @param test_design The test design (optional).
#' @param predictor The predictor used in the multiway classification model (optional).
#' @return A list with class attributes "multiway_classification_result", "classification_result", and "list" containing the observed and predicted values, class probabilities, test design, test indices, and predictor.
#' @inheritParams classification_result
#' @family classification_result
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

 
#' Create a Regression Result Object
#'
#' This function creates a regression result object containing the observed and predicted values, test design, test indices, and predictor.
#'
#' @param observed A vector of observed values.
#' @param predicted A vector of predicted values.
#' @param testind A vector of indices for the test data (optional).
#' @param test_design The test design (optional).
#' @param predictor The predictor used in the regression model (optional).
#' @return A list with class attributes "regression_result", "classification_result", and "list" containing the observed and predicted values, test design, test indices, and predictor.
#' @family classification_result
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

