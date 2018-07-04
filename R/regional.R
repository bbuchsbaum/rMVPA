

combine_regional_results = function(results, ids) {
  if (is.factor(results$result[[1]]$observed)) {
    results %>% dplyr::rowwise() %>% dplyr::do( {

      tib1 <- tibble::tibble(
          ROINUM=rep(.$id, length(.$result$observed)),
          observed=.$result$observed,
          pobserved=sapply(seq_along(.$result$observed), function(i) .$result$probs[i, .$result$observed[i]]),
          predicted=.$result$predicted,
          correct=as.character(.$result$observed) == as.character(.$result$predicted)
      )
      tib2 <- tibble::as_tibble(.$result$probs)
      names(tib2) <- paste0("prob_", names(tib2))
      cbind(tib1, tib2)
    })
  } else {
   
    results %>% dplyr::rowwise() %>% dplyr::do(
      tibble::tibble(
        ROINUM=rep(.$id, length(.$result$observed)),
        observed=.$result$observed,
        predicted=.$result$predicted)
    )
  }
}



#' Region of interest based MVPA analysis
#' 
#' Run a separate MVPA analysis for multiple disjoint regions of interest.
#' 
#' @param model_spec a \code{mvpa_model} instance
#' @param region_mask a \code{BrainVolume} or \code{BrainSurface} where each region is identified by a unique integer. 
#'        Every non-zero set of positive integers will be used to define a set of voxels for clasisifcation analysis.
#' @param return_fits whether to return model fit for every ROI (default is \code{FALSE} to save memory)
#' 
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label 
#'         (e.g. Accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
run_regional <- function(model_spec, region_mask, return_fits=FALSE) {  
  
  ## Get the set of unique ROIs (all unique integers > 0 in provided mask)
  
  region_vec <- as.vector(region_mask)
  region_set <- sort(as.integer(unique(region_vec[region_vec > 0])))
  
  if (length(region_set) < 1) {
    stop("run_regional: invalid ROI mask, number of ROIs = 0")
  }
  
  vox_iter <- lapply(region_set, function(rnum) which(region_vec == rnum))
  lens <- sapply(vox_iter, length)
  
  if (any(lens < 2)) {
    warning(paste("some ROIs have less than two voxels, removing them from list: ", paste(region_set[lens < 2], collapse=",")))
    vox_iter <- vox_iter[lens >= 2]
    region_set <- region_set[lens >= 2]
  }
  
  flog.info("model is: %s", model_spec$model$label)
  
  ## run mvpa for each region
  results <- mvpa_iterate(model_spec, vox_iter, ids=region_set,compute_performance=TRUE, return_fits = return_fits)
  

  ## compile performance results
  perf_mat <- do.call(rbind, results$performance)
  
  ## generate volumetric results
  vols <- lapply(1:ncol(perf_mat), function(i) fill(region_mask, cbind(results$id, perf_mat[,i])))
  names(vols) <- colnames(perf_mat)
  
  perf_mat <- tibble::as_tibble(perf_mat) %>% dplyr::mutate(ROINUM = unlist(results$id)) %>% dplyr::select(ROINUM, dplyr::everything())
  
  ## compile full prediction table
  prediction_table <- combine_regional_results(results, results$id)
  
  fits <- if (return_fits) {
    lapply(results$result, "[[", "predictor")
  }
  
  list(performance_table=perf_mat, prediction_table=prediction_table, vol_results=vols, fits=fits)
}
