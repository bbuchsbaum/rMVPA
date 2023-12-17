

#' Combine Regional Results
#'
#' This function combines regional results from a list into a single data frame.
#'
#' @param results A list of regional results.
#' @return A data frame with combined regional results.
#' @details
#' The function is used to combine regional results from a list into a single data frame.
#' It handles both factor and non-factor observed values and creates a combined data frame with the corresponding columns.
#' @keywords internal
combine_regional_results = function(results) {
  roinum=NULL
  .rownum=NULL
  
  # Check if the observed values are factors (for categorical data)
  if (is.factor(results$result[[1]]$observed)) {
    results %>% dplyr::rowwise() %>% dplyr::do( {
      
      # Create a tibble containing observed, predicted, and additional information
      tib1 <- tibble::tibble(
        .rownum=seq_along(.$result$observed),
        roinum=rep(.$id, length(.$result$observed)),
        observed=.$result$observed,
        pobserved=sapply(seq_along(.$result$observed), function(i) .$result$probs[i, .$result$observed[i]]),
        predicted=.$result$predicted,
        correct=as.character(.$result$observed) == as.character(.$result$predicted)
      )
      
      # Create a tibble with the probabilities for each class
      tib2 <- tibble::as_tibble(.$result$probs, .name_repair=.name_repair)
      names(tib2) <- paste0("prob_", names(tib2))
      
      # Combine tib1 and tib2
      cbind(tib1, tib2)
    })
  } else {
    # For non-factor observed values (for continuous data)
    results %>% dplyr::rowwise() %>% dplyr::do(
      tibble::tibble(
        .rownum=seq_along(.$result$observed),
        roinum=rep(.$id, length(.$result$observed)),
        observed=.$result$observed,
        predicted=.$result$predicted)
    )
  }
}

#' Combine prediction tables
#'
#' Combines multiple prediction tables (e.g., from different models or regions) into a single table.
#' Supports weighted combination and collapsing regions.
#'
#' @param predtabs A list of prediction tables (data frames) to be combined.
#' @param wts A vector of weights, with the same length as \code{predtabs}. Default is equal weights.
#' @param collapse_regions A logical value; if TRUE, regions are collapsed in the final prediction table.
#'
#' @return A combined prediction table (data frame).
#' @import dplyr
#' @importFrom purrr map_df
#' @examples
#' # Create example prediction tables
#' observed = factor(sample(letters[1:2], 10, replace = TRUE))
#' predtab1 <- data.frame(.rownum = 1:10,
#'                        roinum = rep(1, 10),
#'                        observed = observed,
#'                        prob_A = runif(10),
#'                        prob_B = runif(10))
#' predtab2 <- data.frame(.rownum = 1:10,
#'                        roinum = rep(2, 10),
#'                        observed = observed,
#'                        prob_A = runif(10),
#'                        prob_B = runif(10))
#'
#' # Combine the tables
#' combined_table <- combine_prediction_tables(list(predtab1, predtab2))
#' @export
combine_prediction_tables <- function(predtabs, wts=rep(1,length(predtabs)), collapse_regions=FALSE) {
  assert_that(length(predtabs) == length(wts))
  assert_that(sum(wts) > 0)
  assert_that(all(purrr::map_lgl(predtabs, is.data.frame)))
  
  wts <- wts/sum(wts)
  
  .weight <- NULL
  .rownum <- NULL
  roinum <-  NULL
  observed <- NULL
  predicted <- NULL
  
  if (is.character(predtabs[[1]]$observed) || is.factor(predtabs[[1]]$observed)) {
    ## applies constant weight to each table and concatenates
    ptab <- map(seq_along(predtabs), function(i) predtabs[[i]] %>% mutate(.tableid=i, .weight=wts[i])) %>% 
      map_df(bind_rows) %>% as_tibble(.name_repair=.name_repair)
    
    probs <- if (collapse_regions) {
      ptab %>% dplyr::group_by(.rownum,observed) %>% summarise_at(vars(starts_with("prob_")), 
                                                             funs(stats::weighted.mean(., w=.weight)))
    } else {
      ## groups over rownames, condition, and roinum, then compute weighted means of probabilities
      ptab %>% dplyr::group_by(.rownum,observed,roinum) %>% summarise_at(vars(starts_with("prob_")), 
                                                                    funs(stats::weighted.mean(., w=.weight)))
    }
    
    p <- probs %>% ungroup() %>% dplyr::select(dplyr::starts_with("prob_"))
    pmat <- as.matrix(p)
      
    pobserved <- pmat[cbind(seq(1,nrow(probs)),as.integer(probs$observed))]
    mc <- max.col(pmat)
    preds <- levels(probs$observed)[mc]
    
    
    prediction_table <- tibble(
      .rownum=probs$.rownum,
      roinum=if (collapse_regions) 1 else probs$roinum,
      observed=probs$observed,
      pobserved=pobserved,
      predicted=preds,
      correct = predicted == probs$observed,
    ) %>% bind_cols(p)
  } else if (is.numeric(predtabs[[1]]$observed)) {
    stop("combining continuous predictions not implemented")
  }
}


#' Merge regional MVPA results
#'
#' Merge multiple regional MVPA results into a single result.
#'
#' @param x A \code{regional_mvpa_result} object.
#' @param ... Additional \code{regional_mvpa_result} objects to be merged.
#'
#' @return A merged \code{regional_mvpa_result} object.
#' @export
merge_results.regional_mvpa_result <- function(x, ...) {
  rlist <- list(x,...)
  combine_prediction_tables(rlist)
}


#' Create a \code{regional_mvpa_result} instance
#'
#' Constructs a regional MVPA result object that stores the results of MVPA analysis in a specific region.
#'
#' @param model_spec A model specification object.
#' @param performance_table A data frame with performance measures.
#' @param prediction_table A data frame with prediction results.
#' @param vol_results A list of voxel-level results.
#' @param fits Optional model fits.
#'
#' @return A \code{regional_mvpa_result} object.
#' @examples
#' # Create example inputs
#' model_spec <- list(dataset = "Example dataset")
#' performance_table <- data.frame(accuracy = c(0.8, 0.85))
#' prediction_table <- data.frame(observed = factor(rep(letters[1:2], 5)),
#'                                 predicted = factor(rep(letters[1:2], 5)))
#' vol_results <- list(vol1 = "Example vol_result 1", vol2 = "Example vol_result 2")
#' fits <- list(fit1 = "Example fit 1", fit2 = "Example fit 2")
#'
#' # Construct a regional_mvpa_result
#' regional_result <- regional_mvpa_result(model_spec, performance_table,
#'                                         prediction_table, vol_results, fits = fits)
#' @export
regional_mvpa_result <- function(model_spec, performance_table, prediction_table, vol_results, fits=fits) {
  ret <- list(model_spec=model_spec, 
              performance_table=performance_table,
              prediction_table=prediction_table,
              vol_results=vol_results,
              fits=fits)
  
  class(ret) <- c("regional_mvpa_result", "list")
  ret
  
}

#' Prepare regional data for MVPA analysis
#'
#' This function processes the input data and prepares the regions for MVPA analysis by extracting
#' voxel indices for each region of interest (ROI) specified in the region_mask.
#'
#' @param model_spec A model specification object.
#' @param region_mask A mask representing different regions in the brain image.
#'
#' @return A list containing information about the regions for further processing:
#'   * allrois: A vector of unique ROI labels.
#'   * region_vec: A vector representation of the region_mask.
#'   * region_set: A sorted vector of unique ROI labels in the region_mask.
#'   * vox_iter: A list containing voxel indices for each ROI.
#'   * lens: A vector containing the number of voxels in each ROI.
#'   * keep: A logical vector indicating if an ROI should be kept for analysis (those with more than one voxel).
#'
#' @examples
#' # Create example inputs
#' model_spec <- list(dataset = "Example dataset")
#' region_mask <- matrix(c(rep(0, 5), rep(1, 5), rep(2, 5), rep(3, 5)), nrow = 5)
#'
#' # Prepare regional data
#' regional_data <- prep_regional(model_spec, region_mask)
#' @export
prep_regional <- function(model_spec, region_mask) {
  allrois <- sort(unique(region_mask[region_mask>0]))
  region_vec <- as.vector(region_mask)
  region_set <- sort(as.integer(unique(region_vec[region_vec > 0])))
  if (length(region_set) < 1) {
    stop("run_regional: invalid ROI mask, number of ROIs = 0")
  }
  
  vox_iter <- lapply(region_set, function(rnum) which(region_vec == rnum & model_spec$dataset$mask > 0))
  lens <- sapply(vox_iter, length)
  keep <- lens > 1
  
  if (all(!keep)) {
    futile.logger::flog.error("run_regional: no ROIs have more than one voxel.")
    stop()
  }
  
  if (any(lens < 2)) {
    futile.logger::flog.warn(paste("some ROIs have less than two voxels, removing them from list: ", paste(region_set[lens < 2], collapse=",")))
    vox_iter <- vox_iter[keep]
    region_set <- region_set[keep]
  }
  
  list(allrois=allrois, region_vec=region_vec, region_set=region_set,
       vox_iter=vox_iter, lens=lens, keep=keep)
  
}

#' Compile performance results and generate volumetric results
#'
#' This function compiles the performance results from the regional MVPA analysis
#' and generates volumetric results based on the region_mask.
#'
#' @param results A list containing the results from regional MVPA analysis.
#' @param region_mask A mask representing different regions in the brain image.
#'
#' @return A list containing:
#'   * vols: A list of volumetric results for each performance metric.
#'   * perf_mat: A data frame containing the compiled performance results with an added 'roinum' column.
#'
#' @importFrom neuroim2 map_values
#' @keywords internal
comp_perf <- function(results, region_mask) {
  roinum <- NULL
  ## compile performance results
  perf_mat <- do.call(rbind, results$performance)
  
  ## generate volumetric results
  ## TODO fails when region_mask is an logical vol
  vols <- lapply(1:ncol(perf_mat), function(i) map_values(region_mask, 
                                                          cbind(as.integer(results$id), 
                                                                perf_mat[,i])))
  names(vols) <- colnames(perf_mat)
  
  perfmat <- tibble::as_tibble(perf_mat,.name_repair=.name_repair) %>% dplyr::mutate(roinum = unlist(results$id)) %>% dplyr::select(roinum, dplyr::everything())
  list(vols=vols, perf_mat=perfmat)
}

#' Run regional MVPA analysis on a specified MVPA model
#'
#' This function runs a regional MVPA analysis using a specified MVPA model and
#' region mask. The analysis can be customized to return model fits, predictions,
#' and performance measures.
#'
#' @param model_spec An object of type \code{mvpa_model} specifying the MVPA model to be used.
#' @param region_mask A mask representing different regions in the brain image.
#' @param return_fits Whether to return model fit for every ROI (default is \code{FALSE} to save memory).
#' @param return_predictions Whether to return full prediction table with per trial probabilities (can be a large table, set \code{FALSE} to limit memory use).
#' @param compute_performance \code{logical} indicating whether to compute performance measures (e.g. Accuracy, AUC).
#' @param coalesce_design_vars Concatenate additional design variables with output stored in `prediction_table`.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A \code{list} of type \code{regional_mvpa_result} containing a named list of \code{NeuroVol} objects,
#' where each element contains a performance metric and is labeled according to the metric used (e.g. Accuracy, AUC).
#'
#' @import itertools
#' @import foreach
#' @import doParallel
#' @import parallel
#' @export
run_regional.mvpa_model <- function(model_spec, region_mask, return_fits=FALSE, 
                                    return_predictions = TRUE, compute_performance=TRUE, 
                                    coalesce_design_vars=FALSE, ...) {  
  ## to get rid of package check warnings
  roinum=NULL
  ###
  
  prepped <- prep_regional(model_spec, region_mask)
  flog.info("model is: %s", model_spec$model$label)
  
  ## run mvpa for each region
  results <- mvpa_iterate(model_spec, prepped$vox_iter, ids=prepped$region_set,
                          compute_performance=compute_performance, 
                          return_fits = return_fits)

  perf <- if (compute_performance) comp_perf(results, region_mask) else list(vols=list(), perf_mat=tibble())
   
  
  ## compile full prediction table
  prediction_table <- if (return_predictions) {
    combine_regional_results(results) 
  }
  
  if (coalesce_design_vars && return_predictions) {
    prediction_table <- coalesce_join(prediction_table, test_design(model_spec$design), 
                                      by=".rownum")
  }
  
  fits <- if (return_fits) {
    lapply(results$result, "[[", "predictor")
  }
  
  regional_mvpa_result(model_spec=model_spec, performance_table=perf$perf_mat, 
                       prediction_table=prediction_table, vol_results=perf$vols, fits=fits)
}


#' Run regional RSA analysis on a specified RSA model
#'
#' This function runs a regional RSA analysis using a specified RSA model and
#' region mask. The analysis can be customized to return model fits and performance measures.
#'
#' @param model_spec An object of type \code{rsa_model} specifying the RSA model to be used.
#' @param region_mask A mask representing different regions in the brain image.
#' @param return_fits Whether to return model fit for every ROI (default is \code{FALSE} to save memory).
#' @param compute_performance \code{logical} indicating whether to compute performance measures (e.g. Accuracy, AUC).
#' @param regtype The regression method ("pearson", "spearman", "lm", or "rfit").
#' @param distmethod The distance computing method ("pearson" or "spearman").
#' @param coalesce_design_vars Concatenate additional design variables with output stored in `prediction_table`.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A \code{list} of type \code{regional_mvpa_result} containing a named list of \code{NeuroVol} objects,
#' where each element contains a performance metric and is labeled according to the metric used (e.g. Accuracy, AUC).
#'
#' @export
run_regional.rsa_model <- function(model_spec, region_mask, return_fits=FALSE, 
                                   compute_performance=TRUE,regtype=c("pearson", "spearman", "lm", "rfit"), 
                                   distmethod=c("pearson", "spearman"), coalesce_design_vars=FALSE, ...) {  
  regtype <- match.arg(regtype)
  distmethod <- match.arg(distmethod)
  
  prepped <- prep_regional(model_spec, region_mask)
  
  results <- rsa_iterate(model_spec, prepped$vox_iter, ids=prepped$region_set, regtype=regtype, distmethod=distmethod)
  perf <- if (compute_performance) comp_perf(results, region_mask) else list(vols=list(), perf_mat=tibble())
  
  regional_mvpa_result(model_spec=model_spec, performance_table=perf$perf_mat, 
                       prediction_table=NULL, vol_results=perf$vols, fits=NULL)
  
}
  


