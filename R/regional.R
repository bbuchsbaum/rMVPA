

combine_regional_results = function(results, ids) {
  if (is.factor(results$result[[1]]$observed)) {
    results %>% rowwise() %>% do( {
      tib1 <- tibble(
          ROINUM=rep(.$id, length(.$result$observed)),
          observed=.$result$observed,
          predicted=.$result$predicted,
          correct=as.character(.$result$observed) == as.character(.$result$predicted)
      )
      tib2 <- as_tibble(.$result$probs)
      names(tib2) <- paste0("prob_", names(tib2))
      cbind(tib1, tib2)
    })
  } else {
    results %>% rowwise() %>% do( {
      tibble(
        ROINUM=rep(.$id, length(.$result$observed)),
        observed=.$result$observed,
        predicted=.$result$predicted,
      )
    })
  }
}



#' run_regional
#' 
#' Run a separate MVPA analysis for multiple disjoint regions of interest.
#' 
#' @param model_spec a \code{mvpa_model} instance
#' @param region_mask a \code{BrainVolume} where each region is identified by a unique integer. Every non-zero set of positive integers will be used to define a set of voxels for clasisifcation analysis.
#' @param save_predictors whether to return prediction model for every ROI (default is \code{FALSE} to save memory)
#' 
#' @return a named list of \code{BrainVolume} objects, where each name indicates the performance metric and label (e.g. accuracy, AUC)
#' @import itertools 
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
run_regional <- function(model_spec, region_mask, save_predictors=FALSE) {  
  
  ## Get the set of unique ROIs (all unique integers > 0 in provided mask)
  region_set <- sort(as.integer(unique(region_mask[region_mask > 0])))
  
  if (length(region_set) < 1) {
    stop("run_regional: invalid ROI mask, number of ROIs = 0")
  }
  
  vox_iter <- lapply(region_set, function(rnum) which(region_mask == rnum))
  lens <- sapply(vox_iter, length)
  
  if (any(lens < 2)) {
    warning(paste("some ROIs have less than two voxels, removing them from list: ", paste(region_set[lens < 2], collapse=",")))
    vox_iter <- vox_iter[lens >= 2]
  }
  
  ## run mvpa for each region
  results <- mvpa_iterate(model_spec, vox_iter, ids=region_set)
  
  ## compile performance results
  perf_mat <- do.call(rbind, results$performance)
  
  ## generate volumetric results
  vols <- lapply(1:ncol(perf_mat), function(i) fill(region_mask, cbind(results$id, perf_mat[,i])))
  names(vols) <- colnames(perf_mat)
  
  ## compile full prediction table
  prediction_table <- combine_regional_results(results, results$id)
  
  list(performance_table=as_tibble(perf_mat), prediction_table=prediction_table, vol_results=vols)
}
