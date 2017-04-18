


wrap_out <- function(perf_mat, dataset, ids) {
  out <- lapply(1:ncol(perf_mat), function(i)  wrap_output(dataset, perf_mat[,i], ids))
  names(out) <- colnames(perf_mat)
  out
}

#' @importFrom futile.logger flog.error flog.info
#' @importFrom dplyr filter bind_rows
do_randomized <- function(model_spec, radius, niter) {
  ret <- foreach (i = 1:niter) %dopar% {
    flog.info("searchlight iteration: %i", i)
    flog.debug("constructing searchlight.")
    
    slight <- get_searchlight(model_spec$dataset, "randomized", radius)
    flog.debug("done.")
    
    vox_iter <- lapply(slight, function(x) x)
    
    len <- sapply(vox_iter, function(x) attr(x, "length"))
    vox_iter <- vox_iter[len >= 2]
    cind <- sapply(vox_iter, attr, "center.index")
    mvpa_iterate(model_spec, vox_iter, cind)
  }
  
  
  
  results <- dplyr::bind_rows(ret)
  good_results <- results %>% filter(!error)
  bad_results <- results %>% filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    futile.logger::flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    futile.logger::flog.error("no valid results for randomized searchlight, exiting.")
  }
  
  all_ind <- sort(unlist(good_results$indices))
  ind_set <- unique(all_ind)
  ind_count <- table(all_ind)
  
  ncols <- length(good_results$performance[[1]])
  perf_mat <- Matrix::sparseMatrix(i=rep(ind_set, each=ncols), j=rep(1:ncols, length(ind_set)), 
                               x=rep(0, length(ind_set)*ncols), dims=c(length(model_spec$dataset$mask), ncols))
  
  
  for (i in 1:nrow(results)) {
    ind <- good_results$indices[[i]]
    m <- kronecker(matrix(good_results$performance[[i]], 1, ncols), rep(1,length(ind)))
    perf_mat[ind,] <- perf_mat[ind,] + m
  }
  
  perf_mat[ind_set,] <- sweep(perf_mat[ind_set,], 1, as.integer(ind_count), FUN="/")
  colnames(perf_mat) <- names(good_results$performance[[1]])
  wrap_out(perf_mat, model_spec$dataset, ind_set)
    
}


do_standard <- function(model_spec, radius) {
  slight <- get_searchlight(model_spec$dataset, "standard", radius)
  vox_iter <- lapply(slight, function(x) x)
  len <- sapply(vox_iter, function(x) attr(x, "length"))
  vox_iter <- vox_iter[len > 2]
  cind <- sapply(vox_iter, attr, "center.index")
  ret <- mvpa_iterate(model_spec, vox_iter, cind)
  
  good_results <- ret %>% filter(!error)
  bad_results <- ret %>% filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    flog.error("no valid results for randomized searchlight, exiting.")
  }
  
  perf_mat <- good_results %>% dplyr::select(performance) %>% (function(x) do.call(rbind, x[[1]]))
  wrap_out(perf_mat, model_spec$dataset, ret[["id"]])
}


#' run_searchlight
#' 
#' @param model_spec a \code{mvpa_model} instance.
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
run_searchlight <- function(model_spec, radius=8, method=c("randomized", "standard"),  niter=4) {
  
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  if (method == "randomized") {
    assert_that(niter > 1)
  }
  
  flog.info("model is: %s", model_spec$model$label)
  
  res <- if (method == "standard") {
    do_standard(model_spec, radius)    
  } else if (method == "randomized") {
    do_randomized(model_spec, radius, niter)
  } 
  
}