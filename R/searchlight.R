

do_randomized <- function(dataset, model_spec, radius, niter) {
  foreach (i = 1:niter) %dopar% {
    vox_iter <- lapply(neuroim::RandomSearchlight(dataset$mask, radius), function(x) x)
    mvpa_iterate(dataset, vox_iter, model_spec)
    
}


#' run_searchlight
#' 
#' @param dataset a \code{mvpa_dataset} instance.
#' @param model_spec a \code{model_spec} instance.
#' @param radius the searchlight radus in millimeters.
#' @param method the type of searchlight (randomized or standard)
#' @param niter the number of searchlight iterations (used only for 'randomized' method)
#' @return a named list of \code{BrainVolume} objects, where each element contains a performance metric (e.g. AUC or Accuracy) at every voxel location.
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom futile.logger flog.info
#' @export
run_searchlight <- function(dataset, model_spec,radius=8, method=c("randomized", "randomized2", "standard"),  
                             niter=4) {
  stopifnot(niter > 1)
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  flog.info("model is: %s", model$model_name)
  
  res <- if (method == "standard") {
    .doStandard(dataset, model, radius, crossval, class_metrics=class_metrics)    
  } else if (method == "randomized") {
    
    res <- foreach(i = 1:niter) %dopar% {
      flog.info("Running randomized searchlight iteration %s", i)   
      .doRandomized(dataset, model, radius, crossval, class_metrics=class_metrics)
    }
    
    dataset$averageOverIterations(res)
  } else if (method == "randomized2") {
    .doRandomized2(dataset, model, radius, crossval, niter, class_metrics=class_metrics)  
  }
  
}