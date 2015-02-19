#' @export
cleanup <- function(model) {
  UseMethod("cleanup")
}

#' @export
cleanup.default <- function(model) {
  model
}

#' @export
cleanup.mvr <- function(model) {
  model$residuals <- NULL
  model$fitted.values <- NULL
  model$model <- NULL
  model
}

#' @export
cleanup.splsda <- function(model) {
  #model$x <- NULL
  model$y <- NULL
  model
}