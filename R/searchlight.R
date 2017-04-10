

do_randomized <- function(dataset, model_spec, radius, niter) {
  ret <- foreach (i = 1:niter) %dopar% {
    slight <- get_searchlight(dataset, "randomized", radius)
    vox_iter <- lapply(slight, function(x) x)
    cind <- sapply(vox_iter, attr, "center.index")
    mvpa_iterate(dataset, vox_iter, model_spec, cind)
  }
  
  results <- dplyr::bind_rows(ret)
  
  all_ind <- sort(unlist(results$indices))
  ind_set <- unique(all_ind)
  
  ncols <- length(results$performance[[1]])
  omat <- Matrix::sparseMatrix(i=rep(ind_set, each=ncols), j=rep(1:ncols, length(ind_set)), 
                               x=rep(0, length(ind_set)*ncols), dims=c(length(ind_set), ncols))
  
  for (i in 1:nrow(results)) {
    print(i)
    ind <- results$indices[[i]]
    m <- kronecker(matrix(results$performance[[i]], 1, ncols), rep(1,length(ind)))
    omat[ind,] <- omat[ind,] + m
  }
  
  colnames(omat) <- names(results$performance[[1]])
  list(indices=ind_set, result_mat=omat)
    
}


do_standard <- function(dataset, model_spec, radius) {
  slight <- get_searchlight(dataset, "standard", radius)
  vox_iter <- lapply(slight, function(x) x)
  len <- sapply(vox_iter, length)
  cind <- sapply(vox_iter, attr, "center.index")
  ret <- mvpa_iterate(dataset, vox_iter, model_spec, cind)
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
run_searchlight <- function(dataset, model_spec, radius=8, method=c("randomized", "standard"),  
                             niter=4) {
  stopifnot(niter > 1)
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  flog.info("model is: %s", model$model_name)
  
  res <- if (method == "standard") {
    do_standard(dataset, model, radius)    
  } else if (method == "randomized") {
    do_randomized(dataset, model, radius, niter)
  } 
  
}