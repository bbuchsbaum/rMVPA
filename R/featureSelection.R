

#' @export
FeatureSelector <- function(method, cutoff_type, cutoff_value) {
  ret <- list(
              cutoff_type=cutoff_type,
              cutoff_value=cutoff_value)
  class(ret) <- c(method, "list")
  ret
}

#' @export
selectFeatures <- function(obj, X, Y) {
  UseMethod("selectFeatures")
}

## TODO requires that X has spatial structure.
# @export
# @import sda
#selectFeatures.searchlight <- function(obj, X, Y) {
#
#}


#featureMask <- function()

#' @export
#' @import sda
selectFeatures.catscore <- function(obj, X, Y) {
  message("selecting features via catscore")
  sda.1 <- sda.ranking(X, Y)
  
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
#' @importFrom assertthat assert_that
selectFeatures.FTest <- function(obj, X, Y) {
  message("selecting features via FTest")
  message("cutoff type", obj$cutoff_type)
  message("cutoff value", obj$cutoff_value)
  
  assertthat::assert_that(obj$cutoff_type %in% c("topk", "top_k", "topp", "top_p"))
  
  ## TODO speed this up... with matrix operations
  pvals <- unlist(lapply(1:ncol(X), function(i) {
    oneway.test(X[,i] ~ Y)$p.value
  }))
  
  
  keep.idx <- if (obj$cutoff_type == "top_k" || obj$cutoff_type == "topk") {
    k <- min(ncol(X), obj$cutoff_value)
    message("k =", k)
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
 selectFeatures.catscore_FTest <- function(obj, X, Y) {
   message("selecting features via catscore_FTest")
   
   logpvals <- unlist(lapply(1:ncol(X), function(i) {
     -log(oneway.test(X[,i] ~ Y)$p.value)
   }))
   
   sda.1 <- sda.ranking(X, Y)
   idx <- sda.1[,1]
   scores <- numeric(length(idx))
   scores[idx] <- sda.1[,2]
   
   composite <- scale(scores) + scale(logpvals)
   message(cor(scores, logpvals))
 
   keep.idx <- if (obj$cutoff_type == "top_k") {
     k <- min(ncol(X), obj$cutoff_value)
     order(composite, decreasing=TRUE)[1:k]
   } else if (obj$cutoff_type == "top_p") {
     if (obj$cutoff_value <= 0 || obj$cutoff_value > 1) {
       stop("selectFeatures.catscoreFTest: with top_p, cutoff_value must be > 0 and <= 1")
     }
     k <- obj$cutoff_value * ncol(X)
     order(composite, decreasing=TRUE)[1:k]
   } else {
     stop(paste("selectFeatures.catscoreFTest: unsupported cutoff_type: ", obj$cutoff_type))
   }
   
   keep <- logical(ncol(X))
   keep[keep.idx] <- TRUE
   
   message("retaining ", sum(keep), "features in matrix with", ncol(X), "columns")
   keep
 }