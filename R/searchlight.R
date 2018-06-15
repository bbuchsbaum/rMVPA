

#' @keywords internal
wrap_out <- function(perf_mat, dataset, ids) {
  out <- lapply(1:ncol(perf_mat), function(i)  wrap_output(dataset, perf_mat[,i], ids))
  names(out) <- colnames(perf_mat)
  out
}

#' @keywords internal
#' @importFrom futile.logger flog.error flog.info
#' @importFrom dplyr filter bind_rows
do_randomized <- function(model_spec, radius, niter, mvpa_fun=mvpa_iterate, pooled_estimate=TRUE, ...) {
 
  ret <- lapply(1:niter, function(i) {
    flog.info("searchlight iteration: %i", i)
    flog.debug("constructing searchlight.")
    
    slight <- get_searchlight(model_spec$dataset, "randomized", radius)
    
    vox_iter <- lapply(slight, function(x) x)
    
    len <- sapply(vox_iter, function(x) attr(x, "length"))
    vox_iter <- vox_iter[len >= 2]
    cind <- sapply(vox_iter, attr, "center.index")
    mvpa_fun(model_spec, vox_iter, cind,...)
  })
  

  results <- dplyr::bind_rows(ret)
  good_results <- results %>% dplyr::filter(!error)
  bad_results <- results %>% dplyr::filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    futile.logger::flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    futile.logger::flog.error("no valid results for randomized searchlight, exiting.")
  }
  
  all_ind <- sort(unlist(good_results$indices))
  ind_count <- table(all_ind)
  ind_set <- unique(all_ind)
  
  if (pooled_estimate) {
    indmap <- do.call(rbind, lapply(1:nrow(good_results), function(i) {
      ind <- good_results$indices[[i]]
      cbind(i, ind)
    }))
  
    respsets <- split(indmap[,1], indmap[,2])
    
    merged_results <- lapply(respsets, function(r1) {
      first <- r1[1]
      rest <- r1[2:length(r1)]
      z1 <- good_results$result[[first]]
      z2 <- good_results$result[rest]
      ff <- purrr:::partial(merge_results, x=z1)
      do.call(ff, z2)
    })
    
    perf_list <- lapply(merged_results, performance)
    perf_mat <- do.call(rbind, perf_list)
    colnames(perf_mat) <- names(perf_list[[1]])
    wrap_out(perf_mat, model_spec$dataset, ind_set)

  } else {
    ncols <- length(good_results$performance[[1]])
    perf_mat <- Matrix::sparseMatrix(i=rep(ind_set, ncols), j=rep(1:ncols, each=length(ind_set)), 
                               x=rep(0, length(ind_set)*ncols), dims=c(length(model_spec$dataset$mask), ncols))
    for (i in 1:nrow(good_results)) {
      ind <- good_results$indices[[i]]
      m <- kronecker(matrix(good_results$performance[[i]], 1, ncols), rep(1,length(ind)))
      perf_mat[ind,] <- perf_mat[ind,] + m
    }
  
  
    perf_mat[ind_set,] <- sweep(perf_mat[ind_set,], 1, as.integer(ind_count), FUN="/")
    colnames(perf_mat) <- names(good_results$performance[[1]])
    wrap_out(perf_mat, model_spec$dataset, ind_set)
  }
    
}


#' @keywords internal
do_standard <- function(model_spec, radius, mvpa_fun=mvpa_iterate, ...) {
  slight <- get_searchlight(model_spec$dataset, "standard", radius)
  vox_iter <- lapply(slight, function(x) x)
  len <- sapply(vox_iter, function(x) attr(x, "length"))
  vox_iter <- vox_iter[len > 2]
  cind <- sapply(vox_iter, attr, "center.index")
  ret <- mvpa_fun(model_spec, vox_iter, cind,...)

  good_results <- ret %>% dplyr::filter(!error)
  bad_results <- ret %>% dplyr::filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    flog.error("no valid results for randomized searchlight, exiting.")
  }
  
  perf_mat <- good_results %>% dplyr::select(performance) %>% (function(x) do.call(rbind, x[[1]]))
  wrap_out(perf_mat, model_spec$dataset, ret[["id"]])
}


#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom futile.logger flog.info flog.error flog.debug
#' @export
run_searchlight.mvpa_model <- function(model_spec, radius=8, method=c("randomized", "standard"),  niter=4, ...) {
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  if (method == "randomized") {
    assert_that(niter >= 1)
  }
  
  flog.info("model is: %s", model_spec$model$label)
  
  res <- if (method == "standard") {
    flog.info("running standard searchlight with %s radius ", radius)
    do_standard(model_spec, radius)    
  } else if (method == "randomized") {
    flog.info("running randomized searchlight with %s radius and %s iterations", radius, niter)
    do_randomized(model_spec, radius, niter)
  } 
  
}

#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom futile.logger flog.info flog.error flog.debug
#' @export
run_searchlight.rsa_model <- function(model_spec, radius=8, method=c("randomized", "standard"),  niter=4, 
                                      distmethod=c("spearman", "person"), 
                                      regtype=c("lm", "rfit", "pearson", "spearman"),...) {
  
  regtype <- match.arg(regtype)
  distmethod <- match.arg(distmethod)
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  if (method == "randomized") {
    assert_that(niter >= 1)
  }
  
  res <- if (method == "standard") {
    flog.info("running standard RSA searchlight with %s radius ", radius)
    do_standard(model_spec, radius, mvpa_fun=rsa_iterate, regtype, distmethod)    
  } else if (method == "randomized") {
    flog.info("running randomized RSA searchlight with %s radius and %s iterations", radius, niter)
    do_randomized(model_spec, radius, niter, mvpa_fun=rsa_iterate, regtype,distmethod)
  } 
  
}


