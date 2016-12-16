#' @export
load_libs.caret_model_wrapper <- function(x) {
  for (lib in x$model$library) {
    library(lib, character.only = TRUE)
  }
}

#' @export
caret_model_wrapper <- function(model, tune_grid=NULL, custom_performance=NULL) {
  
  if (!is.null(custom_performance)) {
    assert_that(is.function(custom_performance)) 
  }
  
  ret <- list(model=model,
       model_name=model$label,
       tune_grid=tune_grid,
       custom_performance=custom_performance)
  
  class(ret) <- "caret_model_wrapper"
  ret
  
}


mvpa_crossval <- function()
