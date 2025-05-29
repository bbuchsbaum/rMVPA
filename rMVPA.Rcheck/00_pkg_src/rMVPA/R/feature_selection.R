## TODO integrate mlr "filters"


#' Feature Selection Methods
#'
#' @section Methods:
#' Two feature selection methods are available:
#' \describe{
#'   \item{FTest}{One-way ANOVA F-test for each feature}
#'   \item{catscore}{Correlation-adjusted t-scores using sda.ranking}
#' }
#'
#' @section Cutoff Types:
#' Two types of cutoffs are supported:
#' \describe{
#'   \item{top_k/topk}{Select top k features}
#'   \item{top_p/topp}{Select top p percent of features (0 < p <= 1)}
#' }
#'
#' @name feature_selection
NULL


#' @keywords internal
#' @importFrom stats pf
#' @noRd
matrixAnova <- function(Y, x) {
  if (!is.numeric(x)) stop("x must be numeric")
  if (nrow(x) != length(Y)) stop("x and Y must have compatible dimensions")
  if (any(is.na(x)) || any(is.na(Y))) stop("NA values not supported")
  x <- as.matrix(x)
  Y <- as.numeric(Y)
  k <- max(Y)
  ni <- tabulate(Y)
  n <- dim(x)[1]
  sx2 <- colSums(x^2)
  m <- rowsum(x, Y)
  a <- colSums(m^2/ni)
  b <- colSums(m)^2/n
  mst <- (a - b)/(k - 1)
  mse <- (sx2 - a)/(n - k)
  fa <- mst/mse
  pvalue <- pf(fa, k - 1, n - k, lower.tail = FALSE, log.p = FALSE)
  tab <- cbind(fa, pvalue)
  colnames(tab) <- c("Ftest", "pval")
  if (!is.null(colnames(x))) 
    rownames(tab) <- colnames(x)
  tab
  
}



#' Create a feature selection specification
#'
#' This function creates a feature selection specification using the provided
#' method, cutoff type, and cutoff value.
#'
#' @param method The type of feature selection method to use. Supported methods are "FTest" and "catscore".
#' @param cutoff_type The type of threshold used to select features. Supported cutoff types are "top_k" and "top_p".
#' @param cutoff_value The numeric value of the threshold cutoff.
#' @return A list with a class name equal to the \code{method} argument.
#' @details
#' The available feature selection methods are:
#'   - FTest: Computes a one-way ANOVA for every column in the feature matrix.
#'   - catscore: Computes a correlation adjusted t-test for every column in the matrix using \code{sda.ranking} from the \code{sda} package.
#' @examples
#' fsel <- feature_selector("FTest", "top_k", 1000)
#' fsel <- feature_selector("FTest", "top_p", .1)
#' class(fsel) == "FTest"
#' @export
feature_selector <- function(method, cutoff_type, cutoff_value) {
  ret <- list(
              cutoff_type=cutoff_type,
              cutoff_value=cutoff_value)
  class(ret) <- c(method, "feature_selector", "list")
  ret
}

#' @rdname select_features-methods
#' @param ranking.score The feature score to use. Supported scores are "entropy", "avg", or "max". Default is "entropy".
#' @describeIn select_features Feature selection using the CATSCORE method.
#' @export
#' @importFrom sda sda.ranking
select_features.catscore <- function(obj, X, Y,  ranking.score=c("entropy", "avg", "max"),...) {
  assertthat::assert_that(obj$cutoff_type %in% c("topk", "top_k", "topp", "top_p"))
  ranking.score <- match.arg(ranking.score)
  message("selecting features via catscore")
  
  if (is.numeric(Y)) {
    medY <- median(Y)
    Y <- factor(ifelse(Y > medY, "high", "low"))
  }
  
  
  sda.1 <- sda.ranking(as.matrix(X), Y, ranking.score=ranking.score, fdr=FALSE, verbose=FALSE)
  
  keep.idx <- if (obj$cutoff_type == "top_k") {
    k <- min(ncol(X), obj$cutoff_value)
    sda.1[, "idx"][1:k]
  } else if (obj$cutoff_type == "top_p") {
    if (obj$cutoff_value <= 0 || obj$cutoff_value > 1) {
      stop("select_features.catscore: with top_p, cutoff_value must be > 0 and <= 1")
    }
    k <- max(obj$cutoff_value * ncol(X),1)
    sda.1[, "idx"][1:k]
   
  } else {
    stop(paste("select_features.catscore: unsupported cutoff_type: ", obj$cutoff_type))
  }
  
  
  keep <- logical(ncol(X))
  keep[keep.idx] <- TRUE
  message("retaining ", sum(keep), " features in matrix with ", ncol(X), " columns")
  keep
   
}



#' @rdname select_features-methods
#' @describeIn select_features Feature selection using the F-test method.
#' @export
#' @importFrom assertthat assert_that
select_features.FTest <- function(obj, X, Y,...) {
  message("selecting features via FTest")
  message("cutoff type ", obj$cutoff_type)
  message("cutoff value ", obj$cutoff_value)
  
  assertthat::assert_that(obj$cutoff_type %in% c("topk", "top_k", "topp", "top_p"))
  
  if (is.numeric(Y)) {
    medY <- median(Y)
    Y <- factor(ifelse(Y > medY, "high", "low"))
  }
  
  # Ensure X is numeric
  if (!is.numeric(X)) {
    X <- as.matrix(X)
    if (!is.numeric(X)) {
      stop("X must be convertible to a numeric matrix")
    }
  }
  
  pvals <- matrixAnova(Y, X)[,2]
  
  keep.idx <- if (obj$cutoff_type == "top_k" || obj$cutoff_type == "topk") {
    k <- min(ncol(X), obj$cutoff_value)
    order(pvals)[1:k]
  } else if (obj$cutoff_type == "top_p" || obj$cutoff_type == "topp") {
    if (obj$cutoff_value <= 0 || obj$cutoff_value > 1) {
      stop("select_features.FTest: with top_p, cutoff_value must be > 0 and <= 1")
    }
    k <- max(ceiling(obj$cutoff_value * ncol(X)), 1)
    order(pvals)[1:k]
  } else {
  
    stop(paste("select_features.FTest: unsupported cutoff_type: ", obj$cutoff_type))
  }
  
  
  keep <- logical(ncol(X))
  keep[keep.idx] <- TRUE
  
  message("retaining ", sum(keep), " features in matrix with ", ncol(X), " columns")
  
  keep
  
}



# Common validation function
validate_cutoff <- function(type, value, ncol) {
  type <- tolower(type)
  if (!type %in% c("top_k", "topk", "top_p", "topp")) {
    stop("Cutoff type must be one of: top_k, topk, top_p, topp")
  }
  
  if (grepl("p$", type)) {
    if (value <= 0 || value > 1) {
      stop("For percentage cutoff, value must be > 0 and <= 1")
    }
    max(ceiling(value * ncol), 1)
  } else {
    min(value, ncol)
  }
}


#' @export
#' @method print feature_selector
print.feature_selector <- function(x, ...) {
  cat("Feature Selector Object\\n")
  cat("-----------------------\\n")
  cat("Method:        ", class(x)[1], "\\n")
  cat("Cutoff Type:   ", x$cutoff_type, "\\n")
  cat("Cutoff Value:  ", x$cutoff_value, "\\n")
  invisible(x)
}


