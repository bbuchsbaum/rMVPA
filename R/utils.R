
#' groupMeans
#' compute the average vector, splitting a matrix for each level of a grouping variable
#' @param X the matrix
#' @param margin the margin to average over (1 for rows, 2 for columns)
#' @export
groupMeans <- function(X, margin, group) {
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

