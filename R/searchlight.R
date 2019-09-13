
#' @keywords internal
wrap_out <- function(perf_mat, dataset, ids) {
  out <- lapply(1:ncol(perf_mat), function(i)  wrap_output(dataset, perf_mat[,i], ids))
  names(out) <- colnames(perf_mat)
  out
}

#' @keywords internal
combine_randomized <- function(model_spec, good_results, bad_results) {
  all_ind <- sort(unlist(good_results$indices))
  ind_count <- table(all_ind)
  ind_set <- unique(all_ind)
  ncols <- length(good_results$performance[[1]])
  
  perf_mat <- Matrix::sparseMatrix(i=rep(ind_set, ncols), j=rep(1:ncols, each=length(ind_set)), 
                                   x=rep(0, length(ind_set)*ncols), 
                                   dims=c(length(model_spec$dataset$mask), ncols))
  
  for (i in 1:nrow(good_results)) {
    ind <- good_results$indices[[i]]
    m <- kronecker(matrix(good_results$performance[[i]], 1, ncols), rep(1,length(ind)))
    perf_mat[ind,] <- perf_mat[ind,] + m
  }
  
  perf_mat[ind_set,] <- sweep(perf_mat[ind_set,], 1, as.integer(ind_count), FUN="/")
  colnames(perf_mat) <- names(good_results$performance[[1]])
  wrap_out(perf_mat, model_spec$dataset, ind_set)
}

# pool classiifer results collected over a set of overlapping indices
#' @keywords internal
#' @noRd
pool_results <- function(good_results) {
  all_ind <- sort(unlist(good_results$indices))
  ind_count <- table(all_ind)
  ind_set <- unique(all_ind)
  
  indmap <- do.call(rbind, lapply(1:nrow(good_results), function(i) {
    ind <- good_results$indices[[i]]
    cbind(i, ind)
  }))
  
  respsets <- split(indmap[,1], indmap[,2])
  
  merged_results <- furrr::future_map(respsets, function(r1) {
    if (length(r1) > 1) {
      first <- r1[1]
      rest <- r1[2:length(r1)]
      z1 <- good_results$result[[first]]
      z2 <- good_results$result[rest]
      ff <- purrr::partial(merge_results, x=z1)
      do.call(ff, z2)
    } else {
      good_results$result[[r1[1]]]
    }
  })
}

#' @keywords internal
pool_randomized <- function(model_spec, good_results, bad_results) {
  if (nrow(good_results) == 0) {
    stop("searchlight: no searchlight samples produced valid results")
  }
  
  merged_results <- pool_results(good_results)
  perf_list <- furrr::future_map(merged_results, function(res) compute_performance(model_spec, res))
  
  ind_set <- unique(sort(unlist(good_results$indices)))
  #ind_set <- unique(all_ind)
  
  ncols <- length(perf_list[[1]])
  pmat <- do.call(rbind, perf_list)
  
  perf_mat <- Matrix::sparseMatrix(i=rep(ind_set, ncols), j=rep(1:ncols, each=length(ind_set)), 
                                   x=as.vector(pmat), dims=c(length(model_spec$dataset$mask), ncols))
  
  
  colnames(perf_mat) <- names(perf_list[[1]])
  wrap_out(perf_mat, model_spec$dataset, ids=NULL)
}

#' @keywords internal
#' @importFrom futile.logger flog.error flog.info
#' @importFrom dplyr filter bind_rows
#' @importFrom furrr future_map
do_randomized <- function(model_spec, radius, niter, mvpa_fun=mvpa_iterate, combiner=pool_randomized, ...) {
  error=NULL 
  
  ret <- furrr::future_map(1:niter, function(i) {

    futile.logger::flog.info("searchlight iteration: %s", i)
    futile.logger::flog.info("constructing searchlight.")
    
    slight <- get_searchlight(model_spec$dataset, "randomized", radius)
  
    cind <- purrr::map_int(slight, ~ .@parent_index)
    mvpa_fun(model_spec, slight, cind,...)
  })
  
  nmodels <- sum(unlist(sapply(ret, nrow)))
  futile.logger::flog.info("number of models fit: %s", nmodels)
 
  results <- dplyr::bind_rows(ret)
  
  good_results <- results %>% dplyr::filter(error == FALSE)
  bad_results <- results %>% dplyr::filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    futile.logger::flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    futile.logger::flog.error("no valid results for randomized searchlight, exiting.")
  }
  
  ## could simply merge all searchlights to produce global classification measure  
  combiner(model_spec, good_results)
}

#' @keywords performance
combine_standard <- function(model_spec, good_results, bad_results) {
  ind <- unlist(good_results$id)
  perf_mat <- good_results %>% dplyr::select(performance) %>% (function(x) do.call(rbind, x[[1]]))
  wrap_out(perf_mat, model_spec$dataset, ind)
}


#' @keywords internal
do_standard <- function(model_spec, radius, mvpa_fun=mvpa_iterate, combiner=combine_standard, ...) {
  error=NULL
  flog.info("creating standard searchlight")
  slight <- get_searchlight(model_spec$dataset, "standard", radius)
  cind <- which(model_spec$dataset$mask > 0)
  flog.info("running standard searchlight iterator")
  ret <- mvpa_fun(model_spec, slight, cind,...)

  good_results <- ret %>% dplyr::filter(!error)
  bad_results <- ret %>% dplyr::filter(error == TRUE)
  
  if (nrow(bad_results) > 0) {
    flog.info(bad_results$error_message)
  }
  
  if (nrow(good_results) == 0) {
    ## TODO print out some debug information
    flog.error("no valid results for standard searchlight, exiting.")
  }
  
  combiner(model_spec, good_results, bad_results)
}


#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom futile.logger flog.info flog.error flog.debug
#' @param combiner a function that combines results into an appropriate output.
#' 
#' @references 
#' Bjornsdotter, M., Rylander, K., & Wessberg, J. (2011). A Monte Carlo method for locally multivariate brain mapping. Neuroimage, 56(2), 508-516.
#' 
#' Kriegeskorte, N., Goebel, R., & Bandettini, P. (2006). Information-based functional brain mapping. Proceedings of the National academy of Sciences of the United States of America, 103(10), 3863-3868.
#' @export
#' @importFrom purrr pmap
#' @rdname run_searchlight
#' @examples 
#'  
#' dataset <- gen_sample_dataset(c(4,4,4), 100, blocks=3)
#' cval <- blocked_cross_validation(dataset$design$block_var)
#' model <- load_model("sda_notune")
#' mspec <- mvpa_model(model, dataset$dataset, design=dataset$design, model_type="classification", crossval=cval)
#' res <- run_searchlight(mspec,radius=8, method="standard")
#' 
#' # a custom "combiner" can be used to post-process the output of the searchlight classifier for special cases.
#' # in the example below the supplied "combining function" extracts the predicted probability of the correct class 
#' # for every voxel and every trial and then stores them in a data.frame.
#' 
#' \dontrun{ 
#' custom_combiner <- function(mspec, good, bad) { 
#'    good %>% pmap(function(result, id,...) { 
#'      data.frame(trial=1:length(result$observed), id=id, prob=prob_observed(result)) 
#'    }) %>% bind_rows()
#' }
#' 
#' res2 <- run_searchlight(mspec,radius=8, method="standard", combiner=custom_combiner)
#' }
#' 
run_searchlight.mvpa_model <- function(model_spec, radius=8, method=c("randomized", "standard"),  niter=4, combiner=NULL, ...) {
  
  if (radius < 1 || radius > 100) {
    stop(paste("radius", radius, "outside allowable range (1-100)"))
  }
  
  method <- match.arg(method)
  
  if (method == "randomized") {
    assert_that(niter >= 1)
  }
  
  flog.info("model is: %s", model_spec$model$label)
  
  res <- if (method == "standard") {
    if (is.null(combiner)) {
      combiner <- combine_standard
    }
    
    flog.info("running standard searchlight with %s radius ", radius)
    do_standard(model_spec, radius, combiner=combiner)    
  } else if (method == "randomized") {
    
    if (is.null(combiner)) {
      combiner <- pool_randomized
    }
    
    flog.info("running randomized searchlight with %s radius and %s iterations", radius, niter)
    do_randomized(model_spec, radius, niter, combiner=combiner, compute_performance=FALSE)
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
    do_standard(model_spec, radius, mvpa_fun=rsa_iterate, combiner=combine_standard, regtype, distmethod)    
  } else if (method == "randomized") {
    flog.info("running randomized RSA searchlight with %s radius and %s iterations", radius, niter)
    do_randomized(model_spec, radius, niter, mvpa_fun=rsa_iterate, combiner=combine_randomized, regtype,distmethod)
  } 
  
}