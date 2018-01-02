





#' rsa_model
#' 
#' @param Rmat a square relationship matrix with dimensions equal to number of trials in design.
#' @param is_dist \code{TRUE} if the \code{Rmat} a distance matrix
#' @param dataset a \code{mvpa_dataset} instance
#' @param design a \code{rsa_design} instance
#' @param model a estimation model
#' @param crossval a \code{cross_validation} instance
#' @param performance an optional custom function for computing performance metrics.
#' @export
rsa_model <- function(Rmat, 
                      is_dist=TRUE,
                      dataset,
                      design,
                      crossval, 
                      tune_grid=NULL, 
                      performance=NULL,
                      class_metrics=TRUE) {
  
  stop("not implemented")
  
  
}
