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
  # Ensure observed and predicted share the same factor levels so that
  # operations like equality comparisons work without errors.
  observed <- as.factor(observed)
  predicted <- factor(predicted, levels = levels(observed))

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
#' @rdname sub_result-methods
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
#' @rdname sub_result-methods
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

#' @export
#' @method print classification_result
print.classification_result <- function(x, ...) {
  UseMethod("print")
}

#' @export
#' @method print regression_result
print.regression_result <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  stat_style <- crayon::italic$blue
  
  # Print header
  cat("\n", header_style("Regression Result"), "\n\n")
  
  # Basic information
  cat(section_style("- Data Summary"), "\n")
  cat(info_style("  - Observations: "), number_style(length(x$observed)), "\n")
  cat(info_style("  - Test Indices: "), 
      if(!is.null(x$testind)) number_style(length(x$testind)) else crayon::red("None"), "\n")
  
  # Performance metrics
  cat(section_style("- Performance Metrics"), "\n")
  mse <- mean((x$observed - x$predicted)^2)
  rmse <- sqrt(mse)
  r2 <- cor(x$observed, x$predicted)^2
  mae <- mean(abs(x$observed - x$predicted))
  
  cat(info_style("  - MSE: "), number_style(sprintf("%.4f", mse)), "\n")
  cat(info_style("  - RMSE: "), number_style(sprintf("%.4f", rmse)), "\n")
  cat(info_style("  - MAE: "), number_style(sprintf("%.4f", mae)), "\n")
  cat(info_style("  - R2: "), number_style(sprintf("%.4f", r2)), "\n\n")
}

#' @export
#' @method print binary_classification_result
print.binary_classification_result <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  level_style <- crayon::blue
  
  # Print header
  cat("\n", header_style("Binary Classification Result"), "\n\n")
  
  # Basic information
  cat(section_style("- Data Summary"), "\n")
  cat(info_style("  - Observations: "), number_style(length(x$observed)), "\n")
  cat(info_style("  - Classes: "), level_style(paste(levels(x$observed), collapse=", ")), "\n")
  cat(info_style("  - Test Indices: "), 
      if(!is.null(x$testind)) number_style(length(x$testind)) else crayon::red("None"), "\n")
  
  # Performance metrics
  cat(section_style("- Performance Metrics"), "\n")
  conf_mat <- table(Observed=x$observed, Predicted=x$predicted)
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  sensitivity <- conf_mat[2,2] / sum(conf_mat[2,])
  specificity <- conf_mat[1,1] / sum(conf_mat[1,])
  
  cat(info_style("  - Accuracy: "), number_style(sprintf("%.4f", accuracy)), "\n")
  cat(info_style("  - Sensitivity: "), number_style(sprintf("%.4f", sensitivity)), "\n")
  cat(info_style("  - Specificity: "), number_style(sprintf("%.4f", specificity)), "\n\n")
}

#' @export
#' @method print multiway_classification_result
print.multiway_classification_result <- function(x, ...) {
  # Ensure crayon is available
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty printing. Please install it.")
  }
  
  # Define color scheme
  header_style <- crayon::bold$cyan
  section_style <- crayon::yellow
  info_style <- crayon::white
  number_style <- crayon::green
  level_style <- crayon::blue
  
  # Print header
  cat("\n", header_style("Multiway Classification Result"), "\n\n")
  
  # Basic information
  cat(section_style("- Data Summary"), "\n")
  cat(info_style("  - Observations: "), number_style(length(x$observed)), "\n")
  cat(info_style("  - Number of Classes: "), number_style(length(levels(x$observed))), "\n")
  cat(info_style("  - Classes: "), level_style(paste(levels(x$observed), collapse=", ")), "\n")
  cat(info_style("  - Test Indices: "), 
      if(!is.null(x$testind)) number_style(length(x$testind)) else crayon::red("None"), "\n")
  
  # Performance metrics
  cat(section_style("- Performance Metrics"), "\n")
  conf_mat <- table(Observed=x$observed, Predicted=x$predicted)
  accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
  
  # Calculate per-class metrics
  class_metrics <- lapply(levels(x$observed), function(cls) {
    tp <- sum(x$observed == cls & x$predicted == cls)
    total <- sum(x$observed == cls)
    recall <- tp / total
    precision <- tp / sum(x$predicted == cls)
    f1 <- 2 * (precision * recall) / (precision + recall)
    c(recall=recall, precision=precision, f1=f1)
  })
  names(class_metrics) <- levels(x$observed)
  
  cat(info_style("  - Overall Accuracy: "), number_style(sprintf("%.4f", accuracy)), "\n")
  cat(info_style("  - Per-Class Metrics:"), "\n")
  
  for(cls in levels(x$observed)) {
    metrics <- class_metrics[[cls]]
    cat(info_style("    - "), level_style(cls), ":\n")
    cat(info_style("      - Recall: "), number_style(sprintf("%.4f", metrics["recall"])), "\n")
    cat(info_style("      - Precision: "), number_style(sprintf("%.4f", metrics["precision"])), "\n")
    cat(info_style("      - F1: "), number_style(sprintf("%.4f", metrics["f1"])), "\n")
  }
  cat("\n")
}
