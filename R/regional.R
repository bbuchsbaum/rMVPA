#' @keywords internal
get_unique_regions.NeuroVol <- function(region_mask, ...) {
  sort(unique(region_mask[region_mask > 0]))
}

#' @keywords internal
get_unique_regions.NeuroSurface <- function(region_mask, ...) {
  sort(unique(region_mask@data[region_mask@data > 0]))
}

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
#' @noRd
combine_regional_results = function(results) {
  roinum=NULL
  .rownum=NULL
  
  # Check if the observed values are factors (for categorical data)
  if (is.factor(results$result[[1]]$observed)) {
    results %>% dplyr::rowwise() %>% dplyr::do( {
      
      # Use test indices if provided in the classification_result to preserve
      # alignment with the original test design (important when some rows are dropped)
      row_index <- if (!is.null(.$result$testind)) .$result$testind else seq_along(.$result$observed)

      # Create a tibble containing observed, predicted, and additional information
      tib1 <- tibble::tibble(
        .rownum=row_index,
        roinum=rep(.$id, length(row_index)),
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
    results %>% dplyr::rowwise() %>% dplyr::do({
      row_index <- if (!is.null(.$result$testind)) .$result$testind else seq_along(.$result$observed)
      tibble::tibble(
        .rownum=row_index,
        roinum=rep(.$id, length(row_index)),
        observed=.$result$observed,
        predicted=.$result$predicted)
    })
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
                                                                  list(~stats::weighted.mean(., w = .weight)))
    } else {
     
      ## groups over rownames, condition, and roinum, then compute weighted means of probabilities
      ptab %>% dplyr::group_by(.rownum,observed,roinum) %>% summarise_at(vars(starts_with("prob_")), 
                                                                         list(~stats::weighted.mean(., w = .weight)))
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


#'
#' Merge multiple regional MVPA results into a single result.
#'
#' @param obj A \code{regional_mvpa_result} object.
#' @param ... Additional \code{regional_mvpa_result} objects to be merged.
#'
#' @return A merged \code{regional_mvpa_result} object.
#' @rdname merge_results-methods
#' @inheritParams merge_results
#' @method merge_results regional_mvpa_result
#' @export
merge_results.regional_mvpa_result <- function(obj, ...) {
  rlist <- list(obj, ...)
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
#' # Create example data
#' sample_data <- gen_sample_dataset(c(5, 5, 5), nobs = 100, blocks = 4)
#'
#' # Create a simple region mask with 3 ROIs
#' mask_vol <- sample_data$dataset$mask
#' region_mask <- neuroim2::NeuroVol(
#'   sample(1:3, size = sum(mask_vol > 0), replace = TRUE),
#'   space = neuroim2::space(mask_vol),
#'   indices = which(mask_vol > 0)
#' )
#'
#' # Create a basic model spec
#' model_spec <- list(dataset = sample_data$dataset)
#'
#' # Prepare regional data
#' regional_data <- prep_regional(model_spec, region_mask)
#' @export
prep_regional <- function(model_spec, region_mask) {
  allrois <- get_unique_regions(region_mask)
  ##allrois <- sort(unique(region_mask[region_mask>0]))
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
#' @noRd
comp_perf <- function(results, region_mask) {
 
  roinum <- NULL
  ## compile performance results
  perf_mat <- tryCatch({
    do.call(rbind, results$performance)
  }, error = function(e) {
    message("Warning: Error creating performance matrix: ", e$message)
    return(NULL)
  })
  
  # Ensure we keep original names, make unique if duplicates exist
  perf_mat <- as_tibble(perf_mat, .name_repair = "unique")
  
  # Check if perf_mat is NULL or has 0 columns
  if (is.null(perf_mat) || !is.data.frame(perf_mat) || ncol(perf_mat) == 0) {
    message("Warning: Performance matrix is empty or invalid. Returning empty results.")
    return(list(vols=list(), perf_mat=tibble::tibble(roinum=unlist(results$id))))
  }
  
  ## generate volumetric results
  ## TODO fails when region_mask is an logical vol
  vols <- lapply(1:ncol(perf_mat), function(i) map_values(region_mask,
                                                         cbind(as.integer(results$id),
                                                               perf_mat[[i]])))
  names(vols) <- names(perf_mat)
  
  perfmat <- tibble::as_tibble(perf_mat,.name_repair=.name_repair) %>% dplyr::mutate(roinum = unlist(results$id)) %>% dplyr::select(roinum, dplyr::everything())
  list(vols=vols, perf_mat=perfmat)
}

#' @rdname run_regional-methods
#' @param compute_performance Logical indicating whether to compute performance metrics (defaults to \code{model_spec$compute_performance}).
#' @param return_predictions Logical indicating whether to combine a full prediction table (defaults to \code{model_spec$return_predictions}).
#' @param return_fits Logical indicating whether to return the fitted models (defaults to \code{model_spec$return_fits}).
#' @details This function serves as the base implementation for regional analyses, orchestrating data preparation, iteration over regions, performance computation, and result aggregation. Specific `run_regional` methods for different model classes may call this function or provide specialized behavior.
#' @export
run_regional_base <- function(model_spec,
                              region_mask,
                              coalesce_design_vars = FALSE,
                              processor = NULL,
                              verbose = FALSE,
                              compute_performance = model_spec$compute_performance,
                              return_predictions = model_spec$return_predictions,
                              return_fits = model_spec$return_fits,
                        
                              ...) {

  
 
  # 1) Prepare regions
  prepped <- prep_regional(model_spec, region_mask)
  
  # 2) Iterate over regions
  results <- mvpa_iterate(
    model_spec,
    prepped$vox_iter,
    ids = prepped$region_set,
    processor = processor,
    verbose = verbose,
    analysis_type = "regional",
    ...
  )
  
  # 3) Performance computation
  perf <- if (isTRUE(compute_performance)) {
    comp_perf(results, region_mask)
  } else {
    list(vols = list(), perf_mat = tibble::tibble())
  }
  
  # 4) Predictions
  prediction_table <- NULL
  if (isTRUE(return_predictions)) {
    prediction_table <- combine_regional_results(results)
    
    if (coalesce_design_vars && !is.null(prediction_table)) {
      prediction_table <- coalesce_join(
        prediction_table,
        test_design(model_spec$design),
        by = ".rownum"
      )
    }
  }
  
  # 5) Fits
  fits <- NULL
  if (isTRUE(return_fits)) {
    fits <- lapply(results$result, "[[", "predictor")
  }
  
  # 6) Construct and return final result
  regional_mvpa_result(
    model_spec        = model_spec,
    performance_table = perf$perf_mat,
    prediction_table  = prediction_table,
    vol_results       = perf$vols,
    fits             = fits
  )
}


#' Default Method for run_regional
#'
#' @rdname run_regional-methods
#' @details This is the fallback method called when no specialized `run_regional` method is found for the class of `model_spec`. It typically calls `run_regional_base`.
#' @export
run_regional.default <- function(model_spec, region_mask, ...) {
  run_regional_base(model_spec, region_mask, ...)
}


#' Regional MVPA for `mvpa_model` Objects
#'
#' @rdname run_regional-methods
#' @details This method provides the standard regional analysis pipeline for objects of class `mvpa_model` by calling `run_regional_base`.
#' @export
run_regional.mvpa_model <- function(model_spec, region_mask,
                                    coalesce_design_vars = FALSE,
                                    processor = NULL,
                                    verbose = FALSE,
                                    ...) {
  
  run_regional_base(
    model_spec,
    region_mask,
    coalesce_design_vars = coalesce_design_vars,
    processor = processor,
    verbose = verbose,
    ...
  )
}


#' Regional MVPA for `rsa_model` Objects
#'
#' @rdname run_regional-methods
#' @param return_fits Whether to return each region's fitted model (default \code{FALSE}).
#' @param compute_performance \code{logical} indicating whether to compute performance metrics (default \code{TRUE}).
#' @details For `rsa_model` objects, `return_predictions` defaults to `FALSE` as standard RSA typically doesn't produce a prediction table in the same way as classification/regression models.
#' @export
run_regional.rsa_model <- function(model_spec, region_mask,
                                   return_fits = FALSE,
                                   compute_performance = TRUE,
                                   coalesce_design_vars = FALSE,
                                   ...) {
  
  run_regional_base(
    model_spec,
    region_mask,
    coalesce_design_vars  = coalesce_design_vars,
    compute_performance   = compute_performance,
    return_fits           = return_fits,
    return_predictions    = FALSE,  # Override default for RSA
    ...
  )
}

#' Regional MVPA for `vector_rsa_model` Objects
#'
#' @rdname run_regional-methods
#' @param return_fits Logical indicating whether to return the fitted models (default \code{FALSE}).
#' @param compute_performance Logical indicating whether to compute performance metrics (default \code{TRUE}).
#' @details For `vector_rsa_model` objects, `return_predictions` defaults to `FALSE` in `run_regional_base`.
#' If `model_spec$return_predictions` is TRUE, this method will assemble an `observation_scores_table`.
#' @importFrom dplyr bind_rows rename mutate row_number left_join
#' @importFrom tidyr unnest
#' @export
run_regional.vector_rsa_model <- function(model_spec, region_mask,
                                         return_fits = FALSE,
                                         compute_performance = TRUE,
                                         coalesce_design_vars = FALSE, # Usually FALSE for RSA
                                         processor = NULL,
                                         verbose = FALSE,
                                         ...) {
  
  # 1) Prepare regions (using base helper)
  prepped <- prep_regional(model_spec, region_mask)
  
  # 2) Iterate over regions using mvpa_iterate
  # The result from merge_results.vector_rsa_model will contain:
  # - performance: list column with the summary performance matrix
  # - result: list column containing list(rsa_scores=scores_vector) or NULL
  iteration_results <- mvpa_iterate(
    model_spec,
    prepped$vox_iter,
    ids = prepped$region_set,
    processor = processor, # Use default processor unless specified
    verbose = verbose,
    analysis_type = "regional",
    ...
  )
  
  # 3) Performance computation (using base helper)
  # This extracts the 'performance' column from iteration_results
  perf <- if (isTRUE(compute_performance)) {
    comp_perf(iteration_results, region_mask)
  } else {
    list(vols = list(), perf_mat = tibble::tibble())
  }
  
  # 4) Assemble observation scores (if requested)
  prediction_table <- NULL
  if (isTRUE(model_spec$return_predictions) && "result" %in% names(iteration_results)) {
    # Filter out NULL results (where return_predictions was FALSE or errors occurred)
    valid_results <- iteration_results[!sapply(iteration_results$result, is.null), ]
    
    if (nrow(valid_results) > 0) {
      # Create a tibble: roinum | rsa_scores_list
      scores_data <- tibble::tibble(
          roinum = valid_results$id, 
          scores_list = lapply(valid_results$result, function(res) res$rsa_scores)
      )
      
      # Unnest to get a long table: roinum | observation_index | rsa_score
      prediction_table <- scores_data %>%
           mutate(observation_index = map(scores_list, seq_along)) %>% # Add observation index within ROI
           tidyr::unnest(cols = c(scores_list, observation_index)) %>% 
           dplyr::rename(rsa_score = scores_list) # Rename the scores column
           
       # Optionally merge design variables (might need adjustment based on score indices)
       if (coalesce_design_vars) {
            # We need a way to map observation_index back to the original design .rownum
            # This assumes scores are in the same order as the original y_train 
            # (which `second_order_similarity` preserves)
            # Need the original design dataframe 
            orig_design <- model_spec$design$design_table # Assuming it's stored here? Check mvpa_design
            if (!is.null(orig_design)) {
                # Add .rownum based on the original sequence
                # This relies on the assumption that the number of scores matches nrow(orig_design)
                num_obs_in_design <- nrow(orig_design)
                prediction_table <- prediction_table %>%
                   # Need to handle potential mismatch if scores length != num_obs_in_design
                   # For now, assume they match and add .rownum directly
                   dplyr::mutate(.rownum = observation_index) %>%
                   # Perform the join
                   coalesce_join(orig_design, by = ".rownum")
            } else {
                 warning("coalesce_design_vars=TRUE but original design table not found in model_spec$design$design_table")
            }
       }
           
    } else {
         warning("return_predictions=TRUE, but no observation scores were returned from processing.")
    }
  }
  
  # 5) Fits (using base logic - check if applicable for vector_rsa)
  # train_model returns scores, not a fit object, so fits will likely be NULL
  fits <- NULL
  if (isTRUE(return_fits)) {
      # The `result` column now holds scores, not fits. This needs reconsideration.
      # fits <- lapply(iteration_results$result, "[[<some_fit_element>") # This won't work
      warning("`return_fits=TRUE` requested for vector_rsa_model, but this model type does not currently return standard fit objects.")
  }
  
  # 6) Construct and return final result (using base constructor)
  regional_mvpa_result(
    model_spec        = model_spec,
    performance_table = perf$perf_mat,
    prediction_table  = prediction_table, # Add the assembled scores table
    vol_results       = perf$vols,
    fits             = fits
  )
}

#' Regional MVPA for `feature_rsa_model` Objects
#'
#' @rdname run_regional-methods
#' @param coalesce_design_vars If \code{TRUE}, merges design variables into the prediction table.
#' @param processor An optional region processor function.
#' @details This method handles regional analysis for `feature_rsa_model` objects, typically calling `run_regional_base`.
#' @export
#' @S3method run_searchlight feature_rsa_model
# run_regional.feature_rsa_model <- function(model_spec, region_mask,
#                                            coalesce_design_vars = FALSE,
#                                            processor = NULL,
#                                            verbose = FALSE,
#                                            ...) {
# 
#   
#   run_regional_base(
#     model_spec,
#     region_mask,
#     coalesce_design_vars = coalesce_design_vars,
#     processor = processor,
#     verbose = verbose,
#     ...
#   )
# }
