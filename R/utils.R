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