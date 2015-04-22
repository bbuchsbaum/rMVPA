

#' @export
FeatureSelector <- function(method, cutoff.type, cutoff.value) {
  ret <- list(
              cutoff.type=cutoff.type,
              cutoff.value=cutoff.value)
  class(ret) <- c(method, "list")
  ret
}

#' @export
selectFeatures <- function(obj, X, Y) {
  UseMethod("selectFeatures")
}


#featureMask <- function()

#' @export
#' @import sda
selectFeatures.catscore <- function(obj, X, Y) {
  message("selecting features via catscore")
  sda.1 <- sda.ranking(X, Y)
  
  keep.idx <- if (obj$cutoff.type == "top_k") {
    k <- min(ncol(X), obj$cutoff.value)
    sda.1[, "idx"][1:k]
  } else if (obj$cutoff.type == "top_p") {
    if (obj$cutoff.value <= 0 || obj$cutoff.value > 1) {
      stop("selectFeatures.catscore: with top_p, cutoff.value must be > 0 and <= 1")
    }
    k <- obj$cutoff.value * ncol(X)
    sda.1[, "idx"][1:k]
   
  } else {
    stop(paste("selectFeatures.catscore: unsupported cutoff.type: ", obj$cutoff.type))
  }
  
  
  keep <- logical(ncol(X))
  keep[keep.idx] <- TRUE
  message("retaining ", sum(keep), "features in matrix with", ncol(X), "columns")
  keep
   
}

#' @export
selectFeatures.FTest <- function(obj, X, Y) {
  message("selecting features via FTest")
 
  
  pvals <- unlist(lapply(1:ncol(X), function(i) {
    oneway.test(X[,i] ~ Y)$p.value
  }))
  
  
  
 
  keep.idx <- if (obj$cutoff.type == "top_k") {
    k <- min(ncol(X), obj$cutoff.value)
    order(pvals)[1:k]
  } else if (obj$cutoff.type == "top_p") {
    if (obj$cutoff.value <= 0 || obj$cutoff.value > 1) {
      stop("selectFeatures.FTest: with top_p, cutoff.value must be > 0 and <= 1")
    }
    k <- obj$cutoff.value * ncol(X)
    order(pvals)[1:k]
  } else {
    stop(paste("selectFeatures.catscore: unsupported cutoff.type: ", obj$cutoff.type))
  }
  
  keep <- logical(ncol(X))
  keep[keep.idx] <- TRUE
  
  message("retaining ", sum(keep), "features in matrix with", ncol(X), "columns")
  
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
 
   keep.idx <- if (obj$cutoff.type == "top_k") {
     k <- min(ncol(X), obj$cutoff.value)
     order(composite, decreasing=TRUE)[1:k]
   } else if (obj$cutoff.type == "top_p") {
     if (obj$cutoff.value <= 0 || obj$cutoff.value > 1) {
       stop("selectFeatures.catscoreFTest: with top_p, cutoff.value must be > 0 and <= 1")
     }
     k <- obj$cutoff.value * ncol(X)
     order(composite, decreasing=TRUE)[1:k]
   } else {
     stop(paste("selectFeatures.catscoreFTest: unsupported cutoff.type: ", obj$cutoff.type))
   }
   
   keep <- logical(ncol(X))
   keep[keep.idx] <- TRUE
   
   message("retaining ", sum(keep), "features in matrix with", ncol(X), "columns")
   keep
 }