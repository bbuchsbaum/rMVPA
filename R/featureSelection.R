

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

#' @export
#' @import sda
selectFeatures.catscore <- function(obj, X, Y) {
  message("selecting features via catscore")
  sda.1 <- sda.ranking(X, Y)
  
  keep <- if (obj$cutoff.type == "top_k") {
    k <- min(ncol(X), obj$cutoff.value)
    idx <- sda.1[, "idx"][1:k]
    keep <- logical(ncol(X))
    keep[idx] <- TRUE
    keep
  } else if (obj$cutoff.type == "top_p") {
    if (obj$cutoff.value <= 0 || obj$cutoff.value > 1) {
      stop("selectFeatures.catscore: with top_p, cutoff.value must be > 0 and <= 1")
    }
    k <- obj$cutoff.value * ncol(X)
    idx <- sda.1[, "idx"][1:k]
    keep <- logical(ncol(X))
    keep[idx] <- TRUE
    keep
  } else {
    stop(paste("selectFeatures.catscore: unsupported cutoff.type: ", obj$cutoff.type))
  }
  
  keep
   
}