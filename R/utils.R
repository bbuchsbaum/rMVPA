
#' group_means
#' 
#' Given a matrix, compute the average vector for each level of a grouping variable.
#' 
#' @param X the matrix
#' @param margin the margin to average over (1 for rows, 2 for columns)
#' @param group the grouping variable, a \code{factor} or \code{integer} vector
#' @export
group_means <- function(X, margin, group) {
  if (margin == 1) {
    xsum <- rowsum(X, group)
    sweep(xsum, 1, table(group), "/") 
  } else if (margin == 2) {
    xsum <- rowsum(t(X), group)
    t(sweep(xsum, 1, table(group), "/"))
  } else {
    stop("'margin' must be 1 or 2")
  }
}

spearman_cor <- function(x, y=NULL, use="everything") {
  cor(x,y,use, method="spearman")
}

kendall_cor <- function(x, y=NULL, use="everything") {
  cor(x,y,use, method="kendall")
}

zeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) == 0)
}


zeroVarianceColumns2 <- function(M) {
  apply(M, 2, sd, na.rm=TRUE) == 0
}

nonzeroVarianceColumns <- function(M) {
  which(apply(M, 2, sd, na.rm=TRUE) > 0)
}

nonzeroVarianceColumns2 <- function(M) {
  apply(M, 2, sd, na.rm=TRUE) > 0
}

removeZeroVarianceColumns <- function(M) {
  noVariance <- which(apply(M, 2, sd, na.rm=TRUE) == 0)
  if (length(noVariance) > 0) {
    M[, hasVariance, drop=FALSE]
  } else {
    M
  }
}
