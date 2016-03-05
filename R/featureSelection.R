
#' FeatureSelector
#' Creates a feature selection specification
#' @param method the type of feature selection
#' @param cutoff_type the type of threshold used to select features.
#' @param cutoff_value the numeric vale of the threshold cutoff
#' @examples 
#' fsel <- FeatureSelector("FTest", "top_k", 1000)
#' fsel <- FeatureSelector("FTest", "top_p", .1)
#' class(fsel) == "FTest"
#' @return a list with class name equal to the \code{method} arguments
#' @details 
#' 
#' The available feature selection methods are:
#' 
#' Ftest: computes a one-way ANOVA for every column in feature matrix
#' catscore: computes a correlation adjusted t-test for every column in matrix using \code{sda.tanking} from he \code{sda} package.
#' 
#' @export
FeatureSelector <- function(method, cutoff_type, cutoff_value) {
  ret <- list(
              cutoff_type=cutoff_type,
              cutoff_value=cutoff_value)
  class(ret) <- c(method, "list")
  ret
}

#' selectFeatures
#' 
#' Given a \code{FeatureSelection} specification object and a dataset return the set of selected features as a binary vector.
#' @param obj the \code{FeatureSelection} object
#' @param X the training features as a \code{matrix}
#' @param Y the dependent variable as a \code{factor} or \code{numeric} variable.
#' @param vox optional coordinates associated with rows of \code{X} matrix.
#' @param additional arguments
#' @return a \code{logical} vector indicating the columns of \code{X} matrix that were selected.
#' @examples 
#' fsel <- FeatureSelector("FTest", "top_k", 10)
#' X <- matrix(rnorm(100*100), 100, 100)
#' Y <- factor(rep(c("a", "b"), each=50))
#' featureMask <- selectFeatures(fsel, X, Y)
#' sum(featureMask) == 10
#' @export
selectFeatures <- function(obj, X, Y, vox, ...) {
  UseMethod("selectFeatures")
}

## TODO requires that X has spatial structure.
# @export
# @import sda
#selectFeatures.searchlight <- function(obj, X, Y, vox) {
#
#}






#' @export
#' @param ranking.score the feature score: entropy, avg, or max.
#' @rdname selectFeatures
#' @import sda
selectFeatures.catscore <- function(obj, X, Y, vox=NULL, ranking.score=c("entropy", "avg", "max")) {
  ranking.score <- match.arg(ranking.score)
  message("selecting features via catscore")
  
  if (is.numeric(Y)) {
    medY <- median(Y)
    Y <- factor(ifelse(Y > medY, "high", "low"))
  }
  
  sda.1 <- sda.ranking(X, Y, ranking.score=ranking.score)
  
  keep.idx <- if (obj$cutoff_type == "top_k") {
    k <- min(ncol(X), obj$cutoff_value)
    sda.1[, "idx"][1:k]
  } else if (obj$cutoff_type == "top_p") {
    if (obj$cutoff_value <= 0 || obj$cutoff_value > 1) {
      stop("selectFeatures.catscore: with top_p, cutoff_value must be > 0 and <= 1")
    }
    k <- obj$cutoff_value * ncol(X)
    sda.1[, "idx"][1:k]
   
  } else {
    stop(paste("selectFeatures.catscore: unsupported cutoff_type: ", obj$cutoff_type))
  }
  
  
  keep <- logical(ncol(X))
  keep[keep.idx] <- TRUE
  message("retaining ", sum(keep), " features in matrix with ", ncol(X), " columns")
  keep
   
}



#selectFeatures.FTest <- function(obj, X, Y) {
#  Y <- as.factor(Y)
#  
#  rsum <- rowsums(X, Y)
#  nlevs <- table(Y)
#  m.i <- sweep(rsum, 1, nlevs, "/")
#  v.i <- aggregate(X, list(Y), var)[-1]
#}
  
  

#' @export
#' @rdname selectFeatures
#' @importFrom assertthat assert_that
selectFeatures.FTest <- function(obj, X, Y, vox=NULL) {
  message("selecting features via FTest")
  message("cutoff type", obj$cutoff_type)
  message("cutoff value", obj$cutoff_value)
  
  assertthat::assert_that(obj$cutoff_type %in% c("topk", "top_k", "topp", "top_p"))
  
  if (is.numeric(Y)) {
    medY <- median(Y)
    Y <- factor(ifelse(Y > medY, "high", "low"))
  }
  
  ## TODO speed this up... with matrix operations
  pvals <- unlist(lapply(1:ncol(X), function(i) {
    oneway.test(X[,i] ~ Y)$p.value
  }))
  
  
  keep.idx <- if (obj$cutoff_type == "top_k" || obj$cutoff_type == "topk") {
    k <- min(ncol(X), obj$cutoff_value)
    order(pvals)[1:k]
  } else if (obj$cutoff_type == "top_p" || obj$cutoff_type == "topp") {
    if (obj$cutoff_value <= 0 || obj$cutoff_value > 1) {
      stop("selectFeatures.FTest: with top_p, cutoff_value must be > 0 and <= 1")
    }
    k <- obj$cutoff_value * ncol(X)
    order(pvals)[1:k]
  } else {
  
    stop(paste("selectFeatures.FTest: unsupported cutoff_type: ", obj$cutoff_type))
  }
  
  message("length(keep.idx)", length(keep.idx))
  
  keep <- logical(ncol(X))
  keep[keep.idx] <- TRUE
  
  message("retaining ", sum(keep), " features in matrix with ", ncol(X), " columns")
  
  keep
  
}


#' @export
# selectFeatures.catscore_FTest <- function(obj, X, Y, vox=NULL) {
#    message("selecting features via catscore_FTest")
#    
#    if (is.numeric(Y)) {
#      medY <- median(Y)
#      Y <- factor(ifelse(Y > medY, "high", "low"))
#    }
#    
#    logpvals <- unlist(lapply(1:ncol(X), function(i) {
#      -log(oneway.test(X[,i] ~ Y)$p.value)
#    }))
#    
#    sda.1 <- sda.ranking(X, Y)
#    idx <- sda.1[,1]
#    scores <- numeric(length(idx))
#    scores[idx] <- sda.1[,2]
#    
#    composite <- scale(scores) + scale(logpvals)
#    message(cor(scores, logpvals))
#  
#    keep.idx <- if (obj$cutoff_type == "top_k") {
#      k <- min(ncol(X), obj$cutoff_value)
#      order(composite, decreasing=TRUE)[1:k]
#    } else if (obj$cutoff_type == "top_p") {
#      if (obj$cutoff_value <= 0 || obj$cutoff_value > 1) {
#        stop("selectFeatures.catscoreFTest: with top_p, cutoff_value must be > 0 and <= 1")
#      }
#      k <- obj$cutoff_value * ncol(X)
#      order(composite, decreasing=TRUE)[1:k]
#    } else {
#      stop(paste("selectFeatures.catscoreFTest: unsupported cutoff_type: ", obj$cutoff_type))
#    }
#    
#    keep <- logical(ncol(X))
#    keep[keep.idx] <- TRUE
#    
#    message("retaining ", sum(keep), "features in matrix with", ncol(X), "columns")
#    keep
#  }