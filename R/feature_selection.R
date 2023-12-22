## TODO integrate mlr "filters"


#' @keywords iternal
#' @importFrom stats pf
matrixAnova <- function(Y, x) {
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



#' Perform feature selection using the CATSCORE method
#'
#' This function selects features from the input data matrix X using the
#' CATSCORE method and the provided feature selection specification.
#'
#' @param obj The feature selection specification created by \code{feature_selector()}.
#' @param X The input data matrix.
#' @param Y The response variable.
#' @param ranking.score The feature score to use. Supported scores are "entropy", "avg", or "max". Default is "entropy".
#' @return A logical vector indicating which features to retain.
#' @details
#' The CATSCORE method computes a correlation adjusted t-test for every column in the matrix using \code{sda.ranking} from the \code{sda} package.
#' @seealso \code{\link{feature_selector}} for creating a feature selection specification.
#' @export
#' @examples
#' fsel <- feature_selector("catscore", "top_k", 1000)
#' X <- as.data.frame(matrix(rnorm(100 * 10), 100, 10))
#' Y <- rep(letters[1:5], 20)
#' selected_features <- select_features(fsel, X, Y, ranking.score = "entropy")
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


#select_features.FisherKernel <- function(obj, ROI, Y, vox, radius=8) {
#  fres <- matrixAnova(Y,X)
#  search <- Searchlight
#}



#' Perform feature selection using the F-test method
#'
#' This function selects features from the input data matrix X using the
#' F-test method and the provided feature selection specification.
#'
#' @param obj The feature selection specification created by \code{feature_selector()}.
#' @param X The input data matrix.
#' @param Y The response variable.
#' @param ... extra args (not used)
#' @return A logical vector indicating which features to retain.
#' @details
#' The F-test method computes a one-way ANOVA for every column in the feature matrix.
#' @seealso \code{\link{feature_selector}} for creating a feature selection specification.
#' @export
#' @examples
#' fsel <- feature_selector("FTest", "top_k", 1000)
#' X <- as.data.frame(matrix(rnorm(100 * 10), 100, 10))
#' Y <- rep(letters[1:5], 20)
#' selected_features <- select_features(fsel, X, Y)
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
  
  pvals <- matrixAnova(Y,X)[,2]
  
  keep.idx <- if (obj$cutoff_type == "top_k" || obj$cutoff_type == "topk") {
    k <- min(ncol(X), obj$cutoff_value)
    order(pvals)[1:k]
  } else if (obj$cutoff_type == "top_p" || obj$cutoff_type == "topp") {
    if (obj$cutoff_value <= 0 || obj$cutoff_value > 1) {
      stop("select_features.FTest: with top_p, cutoff_value must be > 0 and <= 1")
    }
    k <- obj$cutoff_value * ncol(X)
    order(pvals)[1:k]
  } else {
  
    stop(paste("select_features.FTest: unsupported cutoff_type: ", obj$cutoff_type))
  }
  
  
  keep <- logical(ncol(X))
  keep[keep.idx] <- TRUE
  
  message("retaining ", sum(keep), " features in matrix with ", ncol(X), " columns")
  
  keep
  
}


