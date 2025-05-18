#' Run a Custom Analysis Function Regionally
#'
#' Applies a user-defined function to the data within each specified region
#' of interest (ROI) and returns the results as a tibble.
#'
#' @param dataset An `mvpa_dataset` or `mvpa_surface_dataset` object.
#' @param region_mask A `NeuroVol` or `NeuroSurface` object where each region
#'   is identified by a unique integer greater than 0.
#' @param custom_func A function to apply to each ROI's data. It should
#'   accept two arguments:
#'     \itemize{
#'       \item `roi_data`: A matrix or tibble containing the data
#'             (samples x features) for the current ROI.
#'       \item `roi_info`: A list containing `id` (the region number) and
#'             `indices` (the feature indices for this ROI).
#'     }
#'   The function *must* return a named list or a single-row data frame
#'   (or tibble) containing scalar metric values.
#' @param ... Optional arguments passed to `mvpa_iterate` (e.g., `batch_size`).
#' @param .cores Number of cores to use for parallel processing via the
#'   `future` framework. Defaults to 1 (sequential). Set using
#'   `future::plan()` beforehand for more control.
#' @param .verbose Logical. If `TRUE`, prints progress messages during iteration.
#'   Defaults to `FALSE`.
#'
#' @return A `tibble` where each row corresponds to an ROI. It includes:
#'   \itemize{
#'     \item `id`: The ROI identifier (region number).
#'     \item Columns corresponding to the names returned by `custom_func`.
#'     \item `error`: Logical indicating if an error occurred for this ROI.
#'     \item `error_message`: The error message if an error occurred.
#'   }
#'
#' @details
#' This function provides a simplified interface for applying custom analyses
#' per ROI without needing to define a full `mvpa_model` specification or
#' implement S3 methods. It leverages the parallel processing and iteration
#' capabilities of `rMVPA`.
#'
#' The user-supplied `custom_func` performs the core calculation for each
#' ROI. The framework handles extracting data, iterating over ROIs (potentially
#' in parallel), catching errors from `custom_func`, and formatting the
#' output into a convenient flat table.
#'
#' @examples
#' # Generate sample dataset
#' dset_info <- gen_sample_dataset(D = c(8,8,8), nobs = 50, nlevels = 2)
#' dataset_obj <- dset_info$dataset
#' design_obj <- dset_info$design # Not used by custom_func here, but needed for setup
#'
#' # Create a region mask with 3 ROIs
#' mask_arr <- array(0, dim(dataset_obj$mask))
#' mask_arr[1:4, 1:4, 1:4] <- 1
#' mask_arr[5:8, 1:4, 1:4] <- 2
#' mask_arr[1:4, 5:8, 5:8] <- 3
#' region_mask_vol <- NeuroVol(mask_arr, space(dataset_obj$mask))
#'
#' # Define a custom function: calculate mean and sd for each ROI
#' my_roi_stats <- function(roi_data, roi_info) {
#'   # roi_data is samples x features matrix
#'   # roi_info$id is the region number
#'   # roi_info$indices are the feature indices
#'   mean_signal <- mean(roi_data, na.rm = TRUE)
#'   sd_signal <- sd(roi_data, na.rm = TRUE)
#'   num_features <- ncol(roi_data)
#'   list(
#'     roi_id = roi_info$id, # Can include id if desired, or rely on output table
#'     mean_signal = mean_signal,
#'     sd_signal = sd_signal,
#'     n_features = num_features
#'   )
#' }
#'
#' # Run the custom regional analysis
#' \donttest{
#' # Set up parallel processing (optional)
#' 
#' custom_results <- run_custom_regional(dataset_obj, region_mask_vol, my_roi_stats,
#'                                       .cores = 2, .verbose = TRUE)
#' print(custom_results)
#'
#' # Example with an error in one ROI
#' my_error_func <- function(roi_data, roi_info) {
#'   if (roi_info$id == 2) {
#'     stop("Something went wrong in ROI 2!")
#'   }
#'   list(mean_signal = mean(roi_data))
#' }
#'
#' error_results <- run_custom_regional(dataset_obj, region_mask_vol, my_error_func)
#' print(error_results)
#'
#' # Clean up parallel plan
#' future::plan(future::sequential)
#' }
#' @importFrom dplyr bind_rows select rename all_of
#' @importFrom tidyr unnest_wider
#' @importFrom tibble tibble is_tibble
#' @importFrom neuroim2 values indices space
#' @export
run_custom_regional <- function(dataset, region_mask, custom_func, ...,
                                .cores = 1, .verbose = FALSE) {

  # --- Input Validation ---
  if (!inherits(dataset, c("mvpa_dataset", "mvpa_surface_dataset"))) {
    stop("`dataset` must be an 'mvpa_dataset' or 'mvpa_surface_dataset' object.")
  }
  if (!inherits(region_mask, c("NeuroVol", "NeuroSurface"))) {
     stop("`region_mask` must be a 'NeuroVol' or 'NeuroSurface' object.")
  }
  if (!is.function(custom_func)) {
    stop("`custom_func` must be a function.")
  }
  if (!is.numeric(.cores) || .cores < 1 || round(.cores) != .cores) {
      stop("`.cores` must be a positive integer.")
  }


  # --- Setup Parallel Backend ---
  # Note: User is encouraged to set the plan *before* calling for more control
  if (.cores > 1 && !inherits(future::plan(), c("multicore", "multisession", "cluster"))) {
      if (requireNamespace("future.apply", quietly = TRUE)) {
          message("Setting future plan to 'multisession' with ", .cores, " workers for this function call.")
          old_plan <- future::plan(future::multisession, workers = .cores)
          on.exit(future::plan(old_plan), add = TRUE) # Restore previous plan on exit
      } else {
          warning("Parallel execution requested (cores > 1), but 'future' backend is not multisession/multicore ",
                  "and 'future.apply' is not installed to automatically set it. Running sequentially. ",
                  "Use future::plan() to set backend manually.", call. = FALSE)
      }
  }

  # --- Prepare for Iteration ---
  # Create a minimal dummy model spec - needed for mvpa_iterate internals
  dummy_spec <- list(
      dataset = dataset,
      design = NULL, # Not needed for custom func, but mvpa_iterate expects it
      # Key: we tell iterate we want performance metrics computed,
      # as we hijack the 'performance' slot for our custom metrics
      compute_performance = TRUE,
      return_predictions = FALSE, # Not needed
      return_fits = FALSE         # Not needed
  )
  class(dummy_spec) <- c("custom_internal_model_spec", "model_spec", "list") # Basic class

  # Define the internal processor function
  internal_processor <- function(model_spec, roi, rnum, center_global_id = NA) {
      # Extract necessary info for the custom function
      roi_data <- tryCatch({
          neuroim2::values(roi$train_roi) # Assuming train_roi always exists
      }, error = function(e) { NULL }) # Handle cases where ROI data extraction fails

      roi_indices <- tryCatch({
          neuroim2::indices(roi$train_roi)
      }, error = function(e) { integer(0) })

      if (is.null(roi_data)) {
          # If ROI data extraction failed (e.g., ROI empty after filtering)
           return(tibble::tibble(result = list(NULL), indices = list(roi_indices),
                          performance = list(NULL), id = rnum,
                          error = TRUE, error_message = "Failed to extract ROI data",
                          warning = TRUE, warning_message = "Failed to extract ROI data"))
      }

      roi_info <- list(id = rnum, indices = roi_indices)

      tryCatch({
          # Execute the user's custom function
          perf_result_raw <- custom_func(roi_data, roi_info)

          # Validate and format the result
          if (is.data.frame(perf_result_raw) && nrow(perf_result_raw) == 1) {
              perf_list <- as.list(perf_result_raw)
          } else if (is.list(perf_result_raw) && !is.data.frame(perf_result_raw)) {
              if(!all(sapply(perf_result_raw, function(x) length(x) == 1 && is.atomic(x)))) {
                  stop("custom_func must return a named list of scalars or a single-row data.frame/tibble")
              }
              perf_list <- perf_result_raw
          } else {
              stop("custom_func must return a named list of scalars or a single-row data.frame/tibble")
          }

          # Check for unnamed list elements
          if (is.null(names(perf_list)) || any(names(perf_list) == "")) {
              stop("The list or data frame returned by custom_func must have names for all elements/columns.")
          }

          # Wrap into the tibble structure expected by mvpa_iterate
          tibble::tibble(result = list(NULL), # No model result needed
                         indices = list(roi_info$indices),
                         performance = list(perf_list), # Store metrics here
                         id = rnum, error = FALSE, error_message = "~",
                         warning = FALSE, warning_message = "~")

      }, error = function(e) {
          # Handle errors from custom_func
           tibble::tibble(result = list(NULL), indices = list(roi_indices),
                          performance = list(NULL), id = rnum,
                          error = TRUE, error_message = paste("Error in custom_func:", e$message),
                          warning = TRUE, warning_message = paste("Error in custom_func:", e$message))
      })
  }

  # --- Run Iteration ---
  futile.logger::flog.info("Starting custom regional analysis...")
  prepped <- prep_regional(dummy_spec, region_mask)
  iteration_results <- mvpa_iterate(
    dummy_spec,
    prepped$vox_iter,
    ids = prepped$region_set,
    processor = internal_processor,
    verbose = .verbose,
    analysis_type = "regional",
    ... # Pass other mvpa_iterate args
  )
  futile.logger::flog.info("Custom regional analysis iteration complete.")


  # --- Format Final Output ---
  if (nrow(iteration_results) == 0) {
    warning("No ROIs were successfully processed.")
    return(tibble::tibble(id = integer(), error = logical(), error_message = character()))
  }
  
  # Identify expected columns from the first successful result
  first_success_idx <- which(!iteration_results$error)[1]
  expected_names <- if (!is.na(first_success_idx)) {
      names(iteration_results$performance[[first_success_idx]])
  } else {
       # If all errored, there are no expected metric columns
       character(0) 
  }
  
  # Prepare for unnesting - ensure consistent structure even with errors
  results_to_process <- iteration_results %>%
    dplyr::mutate(performance = lapply(seq_len(nrow(iteration_results)), function(i) {
        p <- .data$performance[[i]]
        err <- .data$error[[i]]
        
        # Create a placeholder list with NA for all expected names
        placeholder <- stats::setNames(as.list(rep(NA, length(expected_names))), expected_names)
        
        if (err || is.null(p)) {
          # If error or NULL result, return the placeholder
          placeholder
        } else {
          # If success, fill the placeholder with actual values found
          common_names <- intersect(names(p), expected_names)
          if (length(common_names) > 0) {
              placeholder[common_names] <- p[common_names]
          }
          placeholder
        }
    }))
  

  # Unnest the performance list-column
  final_table <- tryCatch({
      # Ensure performance is treated as a list column for unnesting
      results_to_process$performance <- as.list(results_to_process$performance)
      
      # Use names_repair to handle potential duplicates (though unlikely now)
      tidyr::unnest_wider(results_to_process, "performance", names_repair = "minimal") %>% 
      # Select explicitly to control order and remove intermediate columns
      dplyr::select(dplyr::all_of(c("id", expected_names, "error", "error_message")))
                     
  }, error = function(e){
       warning(paste("Could not automatically flatten performance metrics:", e$message, 
                     "Returning results in a list column."), call. = FALSE)
       # Ensure fallback also has the correct columns, even if performance is just a list
       fallback_names <- c("id", "error", "error_message")
       if ("performance" %in% names(results_to_process)) {
           fallback_names <- c(fallback_names, "performance")
       }
        # Fallback returns the processed structure before unnesting attempt
        results_to_process %>%
         dplyr::select(dplyr::any_of(fallback_names)) # Use any_of for robustness
  })


  futile.logger::flog.info("Finished formatting custom regional results.")
  return(final_table)
}

# Add a dummy method for the internal class to satisfy mvpa_iterate checks
#' @export
#' @keywords internal
process_roi.custom_internal_model_spec <- function(mod_spec, roi, rnum, ...) {
  # This should not be called directly if the processor is provided,
  # but needs to exist.
  stop("Internal error: process_roi called for custom_internal_model_spec")
}


#' Run a Custom Analysis Function in a Searchlight
#'
#' Applies a user-defined function to the data within each searchlight sphere
#' and returns the results, typically as `NeuroVol` or `NeuroSurface` objects
#' within a `searchlight_result` structure.
#'
#' @param dataset An `mvpa_dataset` or `mvpa_surface_dataset` object.
#' @param custom_func A function to apply within each searchlight sphere. It
#'   should accept two arguments:
#'     \itemize{
#'       \item `sl_data`: A matrix or tibble containing the data
#'             (samples x features_in_sphere) for the current sphere.
#'       \item `sl_info`: A list containing information about the sphere,
#'             including `center_index` (the index of the center voxel/vertex),
#'             `indices` (the indices of all features within the sphere), and
#'             potentially `coords` (coordinates of the center).
#'     }
#'   The function *must* return a named list or a single-row data frame
#'   (or tibble) containing scalar metric values. All spheres must return the
#'   same named metrics.
#' @param radius The radius of the searchlight sphere (in mm for volumes,
#'   or vertex connections for surfaces - see `neuroim2::spherical_roi`).
#' @param method The type of searchlight: "standard" (systematically covers
#'   all center voxels) or "randomized" (samples spheres randomly, useful for
#'   large datasets). Defaults to "standard".
#' @param niter The number of iterations for a "randomized" searchlight.
#'   Ignored if `method = "standard"`. Defaults to 100.
#' @param ... Optional arguments passed to `mvpa_iterate` (e.g., `batch_size`).
#' @param .cores Number of cores to use for parallel processing via the
#'   `future` framework. Defaults to 1 (sequential). Set using
#'   `future::plan()` beforehand for more control.
#' @param .verbose Logical. If `TRUE`, prints progress messages during iteration.
#'   Defaults to `FALSE`.
#'
#' @return A `searchlight_result` object (see `rMVPA::wrap_out`). This is a list
#'   containing:
#'   \itemize{
#'     \item `results`: A named list where each element corresponds to a metric
#'           returned by `custom_func`. Each element is itself a
#'           `searchlight_performance` object containing a `NeuroVol` or
#'           `NeuroSurface` (`$data`) with the metric values mapped back to the
#'           brain space, along with summary statistics (`$summary_stats`).
#'     \item `metrics`: A character vector of the metric names.
#'     \item `n_voxels`: total voxels/vertices defined by the mask.
#'     \item `active_voxels`: number of voxels/vertices with results.
#'   }
#'   If `method = "randomized"`, the values in the output maps represent the
#'   average metric value for each voxel across all spheres it participated in.
#'
#' @details
#' This function provides a flexible way to perform custom analyses across the
#' brain using a searchlight approach, without defining a full `mvpa_model`.
#' It handles iterating over searchlight spheres, extracting data, running the
#' custom function (potentially in parallel), handling errors, and combining
#' the results back into brain-space maps.
#'
#' The `custom_func` performs the core calculation for each sphere. The framework
#' manages the iteration, data handling, parallelization, error catching, and
#' result aggregation.
#'
#' For `method = "standard"`, the function iterates through every active voxel/vertex
#' in the dataset mask as a potential sphere center.
#' For `method = "randomized"`, it randomly selects sphere centers for `niter`
#' iterations. The final map represents an average of the results from spheres
#' covering each voxel. This requires the custom function's results to be meaningfully
#' averageable.
#'
#' **Important**: The `custom_func` must consistently return the same set of named
#' scalar metrics for every sphere it successfully processes.
#'
#' @examples
#' # Generate sample dataset
#' dset_info <- gen_sample_dataset(D = c(10, 10, 10), nobs = 30, nlevels = 2)
#' dataset_obj <- dset_info$dataset
#'
#' # Define a custom function: calculate mean and sd within the sphere
#' my_sl_stats <- function(sl_data, sl_info) {
#'   # sl_data is samples x features_in_sphere matrix
#'   # sl_info contains center_index, indices, etc.
#'   mean_signal <- mean(sl_data, na.rm = TRUE)
#'   sd_signal <- sd(sl_data, na.rm = TRUE)
#'   n_features <- ncol(sl_data)
#'   list(
#'     mean_signal = mean_signal,
#'     sd_signal = sd_signal,
#'     n_vox_in_sphere = n_features
#'   )
#' }
#'
#' # Run the custom searchlight (standard method)
#' \donttest{
#' 
#' custom_sl_results <- run_custom_searchlight(dataset_obj, my_sl_stats,
#'                                             radius = 7, method = "standard",
#'                                             .cores = 2, .verbose = TRUE)
#' print(custom_sl_results)
#'
#' # Access the NeuroVol for a specific metric
#' mean_signal_map <- custom_sl_results$results$mean_signal$data
#' # plot(mean_signal_map) # Requires neuroim2 plotting capabilities
#'
#' # Example with an error in some spheres (e.g., if too few voxels)
#' my_error_sl_func <- function(sl_data, sl_info) {
#'   if (ncol(sl_data) < 5) {
#'     stop("Too few voxels in this sphere!")
#'   }
#'   list(mean_signal = mean(sl_data))
#' }
#'
#' error_sl_results <- run_custom_searchlight(dataset_obj, my_error_sl_func,
#'                                            radius = 4, method = "standard")
#' print(error_sl_results) # Errors will be caught, corresponding voxels may be NA
#'
#' # Run randomized searchlight (faster for large datasets/radii)
#' custom_sl_rand_results <- run_custom_searchlight(dataset_obj, my_sl_stats,
#'                                                  radius = 7, method = "randomized",
#'                                                  niter = 50, # Fewer iterations for example
#'                                                  .cores = 2, .verbose = TRUE)
#' print(custom_sl_rand_results)
#'
#' # Clean up parallel plan
#' future::plan(future::sequential)
#' }
#'
#' @importFrom dplyr bind_rows select filter pull mutate
#' @importFrom tidyr unnest_wider
#' @importFrom tibble tibble is_tibble
#' @importFrom Matrix sparseMatrix 
#' @importFrom purrr map map_int map_dbl map_lgl list_assign
#' @importFrom stats setNames sd
#' @importFrom futile.logger flog.info flog.warn flog.error flog.debug
#' @importFrom methods is
#' @importFrom future plan multisession sequential %<-% %globals%
#' @importFrom future.apply future_lapply future_mapply
#' @importFrom assertthat assert_that
#'
#' @seealso \code{\link{run_custom_regional}}, \code{\link{run_searchlight_base}}, \code{\link{get_searchlight}}, \code{\link{mvpa_iterate}}
#' @export
run_custom_searchlight <- function(dataset, custom_func, radius,
                                   method = c("standard", "randomized"),
                                   niter = 100, ...,
                                   .cores = 1, .verbose = FALSE) {

  # --- Input Validation ---
  method <- match.arg(method)
  if (!inherits(dataset, c("mvpa_dataset", "mvpa_surface_dataset"))) {
    stop("`dataset` must be an 'mvpa_dataset' or 'mvpa_surface_dataset' object.")
  }
  if (!is.function(custom_func)) {
    stop("`custom_func` must be a function.")
  }
  if (!is.numeric(radius) || radius <= 0) {
      stop("`radius` must be a positive number.")
  }
  if (method == "randomized" && (!is.numeric(niter) || niter < 1 || round(niter) != niter)) {
      stop("`niter` must be a positive integer for randomized searchlight.")
  }
  if (!is.numeric(.cores) || .cores < 1 || round(.cores) != .cores) {
      stop("`.cores` must be a positive integer.")
  }

  # --- Setup Parallel Backend ---
  if (.cores > 1 && !inherits(future::plan(), c("multicore", "multisession", "cluster"))) {
      if (requireNamespace("future.apply", quietly = TRUE)) {
          message("Setting future plan to 'multisession' with ", .cores, " workers for this function call.")
          old_plan <- future::plan(future::multisession, workers = .cores)
          on.exit(future::plan(old_plan), add = TRUE) # Restore previous plan on exit
      } else {
          warning("Parallel execution requested (cores > 1), but 'future' backend is not multisession/multicore ",
                  "and 'future.apply' is not installed to automatically set it. Running sequentially. ",
                  "Use future::plan() to set backend manually.", call. = FALSE)
      }
  }

  # --- Prepare for Iteration ---
  # Create a minimal dummy model spec - needed for mvpa_iterate internals,
  # but most fields won't be used directly by our custom processor/combiner.
  # Pass dataset for combiner access.
  dummy_spec <- list(
      dataset = dataset,
      design = NULL,
      compute_performance = TRUE, # Hijack performance slot for custom metrics
      return_predictions = FALSE,
      return_fits = FALSE
  )
  class(dummy_spec) <- c("custom_internal_model_spec", "model_spec", "list")

  # Define the internal processor function for ONE searchlight sphere
  internal_processor <- function(model_spec, roi, rnum, center_global_id = NA) {
      # `roi` here is the searchlight object (e.g., ROIVolume)
      # `rnum` is the index of the center voxel/vertex
      
      sl_data <- tryCatch({
          neuroim2::values(roi$train_roi) # Data within the sphere (samples x features)
      }, error = function(e) { NULL })

      sl_indices <- tryCatch({
           neuroim2::indices(roi$train_roi) # Feature indices within the sphere
      }, error = function(e) { integer(0) })

      if (is.null(sl_data) || ncol(sl_data) == 0) {
          # Handle cases where ROI data extraction fails or sphere is empty
           return(tibble::tibble(result = list(NULL), # No model result
                          indices = list(sl_indices), # Indices within sphere
                          performance = list(NULL), # custom_func output goes here
                          id = rnum, # Center voxel index
                          error = TRUE, error_message = "Failed to extract data or sphere empty",
                          warning = TRUE, warning_message = "Failed to extract data or sphere empty"))
      }

      sl_info <- list(center_index = rnum, indices = sl_indices)
      # Could add coords: sl_info$coords <- neuroim2::coords(neuroim2::spatial(roi$train_roi))

      tryCatch({
          # Execute the user's custom function
          perf_result_raw <- custom_func(sl_data, sl_info)

          # Validate and format the result
          if (tibble::is_tibble(perf_result_raw) && nrow(perf_result_raw) == 1) {
              perf_list <- as.list(perf_result_raw)
          } else if (is.list(perf_result_raw) && !is.data.frame(perf_result_raw)) {
              if(!all(sapply(perf_result_raw, function(x) length(x) == 1 && is.atomic(x)))) {
                  stop("custom_func must return a named list of scalars or a single-row data.frame/tibble")
              }
              perf_list <- perf_result_raw
          } else {
              stop("custom_func must return a named list of scalars or a single-row data.frame/tibble")
          }

          # Check for unnamed list elements
          if (is.null(names(perf_list)) || any(names(perf_list) == "")) {
              stop("The list or data frame returned by custom_func must have names for all elements/columns.")
          }
           # Check all are scalar
           if (!all(sapply(perf_list, function(x) length(x) == 1 && is.atomic(x)))) {
              stop("custom_func must return a named list of scalars or a single-row data.frame/tibble (all elements/columns must be single atomic values).")
           }


          # Wrap into the tibble structure expected by mvpa_iterate/combiner
          tibble::tibble(result = list(NULL), # No model result needed
                         indices = list(sl_info$indices), # Indices within sphere
                         performance = list(perf_list), # Store custom metrics here
                         id = rnum, # Center voxel index
                         error = FALSE, error_message = "~",
                         warning = FALSE, warning_message = "~")

      }, error = function(e) {
          # Handle errors from custom_func
          # Ensure the structure is consistent even on error
           tibble::tibble(result = list(NULL), indices = list(sl_indices),
                          performance = list(NULL), # Performance is NULL on error
                          id = rnum,
                          error = TRUE, error_message = paste("Error in custom_func:", e$message),
                          warning = TRUE, warning_message = paste("Error in custom_func:", e$message))
      })
  }

  # --- Run Iteration ---
  flog.info("Starting custom searchlight analysis (method: %s, radius: %s mm)...", method, radius)

  if (method == "standard") {
      slight <- get_searchlight(dataset, "standard", radius)
      center_indices <- which(dataset$mask > 0) # Center on all active voxels
      flog.info("Preparing %d standard searchlight spheres...", length(center_indices))

      iteration_results <- mvpa_iterate(
          dummy_spec,
          slight,
          ids = center_indices,
          processor = internal_processor,
          verbose = .verbose,
          ... # Pass other mvpa_iterate args like batch_size
      )
      
      flog.info("Combining results from standard searchlight...")
      # Pass dataset directly for combiner access
      final_result <- combine_custom_standard(dataset, iteration_results) 

  } else { # method == "randomized"
      # Randomized searchlight needs a loop and a different combiner
       flog.info("Running %d randomized searchlight iterations...", niter)
       all_iteration_results <- list()
       
       for (i in 1:niter) {
            if (.verbose) flog.info("Randomized iteration %d/%d", i, niter)
            # Get random spheres for this iteration
            slight <- get_searchlight(dataset, "randomized", radius) 
            
            # Extract center indices (handling both NeuroVol and NeuroSurface cases)
             center_indices <- if (methods::is(slight[[1]], "ROIVolume")) {
               sapply(slight, function(x) x@parent_index)
             } else if (methods::is(slight[[1]], "ROISurface")) {
               sapply(slight, function(x) x@center_index)
             } else {
                # Fallback/error - should identify centers based on slight type
                 warning("Could not determine center indices for randomized searchlight iteration.")
                 integer(0) 
             }
            
             if (length(slight) == 0 || length(center_indices) == 0) {
                 flog.warn("No searchlight spheres generated in iteration %d. Skipping.", i)
                 next
             }

             iter_res <- mvpa_iterate(
                 dummy_spec,
                 slight,
                 ids = center_indices, # IDs are the center indices for this random batch
                 processor = internal_processor,
                 verbose = FALSE, # Usually too noisy for randomized inner loop
                 ...
             )
             # Store results, including sphere indices for averaging
             # Need indices *within* the sphere, which processor calculates
              all_iteration_results[[length(all_iteration_results) + 1]] <- iter_res
       }
       
       if (length(all_iteration_results) == 0) {
           stop("No results generated from any randomized searchlight iteration.")
       }
       
       # Combine results from all iterations
       combined_iterations <- dplyr::bind_rows(all_iteration_results)
       flog.info("Combining results from randomized searchlight (%d total spheres processed)...", nrow(combined_iterations))
       final_result <- combine_custom_randomized(dataset, combined_iterations)
  }


  flog.info("Finished custom searchlight analysis.")
  return(final_result)
}


#' Combine Custom Standard Searchlight Results
#'
#' Internal function to combine results from a standard custom searchlight run.
#' Creates a `searchlight_result` object with NeuroVol/NeuroSurface for each metric.
#'
#' @param dataset The original mvpa_dataset object.
#' @param iteration_results The raw tibble output from `mvpa_iterate`.
#' @return A `searchlight_result` object.
#' @keywords internal
combine_custom_standard <- function(dataset, iteration_results) {
  good_results <- iteration_results %>% dplyr::filter(!.data$error)
  bad_results <- iteration_results %>% dplyr::filter(.data$error)

  if (nrow(good_results) == 0) {
    flog.error("No successful results for standard custom searchlight. Examining errors:")
    if (nrow(bad_results) > 0) {
      error_summary <- table(bad_results$error_message)
      for (i in seq_along(error_summary)) {
        flog.error("  - %s: %d occurrences", names(error_summary)[i], error_summary[i])
      }
    }
    stop("No valid results for standard custom searchlight: all spheres failed.")
  }

  # Extract performance lists and center IDs
  perf_lists <- good_results$performance
  center_ids <- good_results$id # Center voxel indices for successful spheres

  # Check consistency of metric names across results
  metric_names <- names(perf_lists[[1]])
  if (is.null(metric_names) || any(metric_names == "")) {
       stop("Internal error: First successful result has unnamed/empty metrics.")
  }
  
  # Validate that all successful results have the same metric names
  all_names_consistent <- all(sapply(perf_lists[-1], function(p) {
      identical(names(p), metric_names)
  }))
  if (!all_names_consistent) {
      stop("Custom function returned inconsistent metric names across different spheres.")
  }
  
  num_metrics <- length(metric_names)
  num_results <- nrow(good_results)

  # Create a matrix: rows = successful center voxels, cols = metrics
  perf_mat <- matrix(NA_real_, nrow = num_results, ncol = num_metrics,
                     dimnames = list(NULL, metric_names))

  # Fill the matrix
  for (i in 1:num_results) {
     # Ensure the list has the expected names before assigning
     p_list <- perf_lists[[i]]
     if (identical(names(p_list), metric_names)) {
        perf_mat[i, ] <- unlist(p_list)
     } else {
         # This shouldn't happen due to the check above, but as a safeguard:
         flog.warn("Metric name mismatch for center voxel %d. Setting to NA.", center_ids[i])
         # Attempt partial matching if names are just reordered (less robust)
         # matched_indices <- match(metric_names, names(p_list))
         # perf_mat[i, !is.na(matched_indices)] <- unlist(p_list[matched_indices[!is.na(matched_indices)]])
     }
  }

  # Use the wrap_out structure (adapted)
  out_list <- lapply(1:num_metrics, function(i) {
    metric_name <- metric_names[i]
    metric_data <- perf_mat[, i]
    # Create the performance object, passing the vector and center IDs
    create_searchlight_performance(dataset, metric_data, center_ids) # Needs vector + IDs
  })
  names(out_list) <- metric_names

  # Create the final searchlight_result object
  structure(
    list(
      results = out_list,
      n_voxels = length(dataset$mask),
      active_voxels = sum(dataset$mask > 0),
      metrics = metric_names
    ),
    class = c("searchlight_result", "list")
  )
}


#' Combine Custom Randomized Searchlight Results
#'
#' Internal function to combine results from a randomized custom searchlight run.
#' Averages results for each metric across overlapping spheres.
#'
#' @param dataset The original mvpa_dataset object.
#' @param iteration_results The raw tibble output from *all* iterations of `mvpa_iterate`.
#' @return A `searchlight_result` object.
#' @keywords internal
combine_custom_randomized <- function(dataset, iteration_results) {
   good_results <- iteration_results %>% dplyr::filter(!.data$error)
   bad_results <- iteration_results %>% dplyr::filter(.data$error)

   if (nrow(good_results) == 0) {
      flog.error("No successful results for randomized custom searchlight. Examining errors:")
      if (nrow(bad_results) > 0) {
          error_summary <- table(bad_results$error_message)
          for (i in seq_along(error_summary)) {
              flog.error("  - %s: %d occurrences", names(error_summary)[i], error_summary[i])
          }
      }
      stop("No valid results for randomized custom searchlight: all spheres failed.")
   }

   # Extract performance, sphere indices, and check metric consistency
   perf_lists <- good_results$performance
   sphere_indices_list <- good_results$indices # List of indices *within* each sphere
   
   metric_names <- names(perf_lists[[1]])
    if (is.null(metric_names) || any(metric_names == "")) {
       stop("Internal error: First successful result has unnamed/empty metrics.")
   }
   num_metrics <- length(metric_names)
   
   all_names_consistent <- all(sapply(perf_lists[-1], function(p) {
      identical(names(p), metric_names)
   }))
   if (!all_names_consistent) {
      stop("Custom function returned inconsistent metric names across different spheres.")
   }

   # --- Accumulation Logic (similar to combine_randomized) ---
   
   # Get all unique voxel indices covered by any successful sphere
   all_covered_indices <- unique(sort(unlist(sphere_indices_list)))
   if (length(all_covered_indices) == 0) {
       stop("No voxels were covered by any successful searchlight sphere.")
   }
   
   # Count how many spheres covered each voxel index
   index_counts <- table(unlist(sphere_indices_list))
   # Ensure index_counts uses character keys for subsetting sparseMatrix later
   names(index_counts) <- as.character(names(index_counts)) 
   
   # Create sparse matrices to accumulate results for *each* metric
   # Rows = all voxels in mask, Cols = 1 (for each metric)
   accumulators <- list()
   total_voxels_in_mask <- length(dataset$mask) # Total size needed for sparse matrix

   for (m_name in metric_names) {
       accumulators[[m_name]] <- Matrix::sparseMatrix(
           i = integer(0), j = integer(0), x = numeric(0),
           dims = c(total_voxels_in_mask, 1) 
       )
   }

   # Iterate through each successful sphere result and add its metric values
   # to the corresponding voxels in the accumulators
   for (i in 1:nrow(good_results)) {
       sphere_vox_indices <- sphere_indices_list[[i]]
       perf_values <- perf_lists[[i]] # Named list of metrics for this sphere

       if (length(sphere_vox_indices) > 0 && !is.null(perf_values)) {
           # Add each metric value to the appropriate accumulator
            for (m_name in metric_names) {
                metric_value <- perf_values[[m_name]]
                if (!is.null(metric_value) && is.numeric(metric_value) && length(metric_value) == 1) {
                   # Add the single metric value to all voxels in this sphere
                   tryCatch({
                       # Ensure indices are valid before subsetting
                       valid_indices <- sphere_vox_indices[sphere_vox_indices <= total_voxels_in_mask & sphere_vox_indices > 0]
                       if (length(valid_indices) > 0) {
                          accumulators[[m_name]][valid_indices, 1] <- accumulators[[m_name]][valid_indices, 1] + metric_value
                       }
                    }, error = function(e) {
                         flog.warn("Error adding metric '%s' for sphere %d (center %d): %s", 
                                   m_name, i, good_results$id[i], e$message)
                    })
                } else {
                     flog.warn("Invalid metric value ('%s') for sphere %d (center %d). Skipping.", 
                               m_name, i, good_results$id[i])
                }
            }
       }
   }

   # --- Normalization and Output Wrapping ---
   out_list <- list()
   active_mask_indices_char <- as.character(which(dataset$mask > 0))
   
   for (m_name in metric_names) {
       acc_matrix <- accumulators[[m_name]]
       
       # Identify indices present in both accumulator and counts
       indices_to_normalize_char <- intersect(rownames(acc_matrix), names(index_counts))
       indices_to_normalize_num <- as.numeric(indices_to_normalize_char)

       if (length(indices_to_normalize_num) > 0) {
           # Extract counts for normalization
           counts_for_norm <- index_counts[indices_to_normalize_char]
           
           # Perform normalization (division) using sweep or direct division
           # Ensure counts are numeric for division
            normalized_values <- acc_matrix[indices_to_normalize_num, 1] / as.numeric(counts_for_norm)
           
           # Update the accumulator matrix with normalized values
           acc_matrix[indices_to_normalize_num, 1] <- normalized_values
            
            # Set other values (those with counts but maybe no accumulation?) to 0 or NA? Let's keep them as they are (likely 0)
       } else {
           flog.warn("No overlapping indices found for normalization for metric '%s'.", m_name)
       }

       # Wrap the normalized result vector (column 1 of sparse matrix)
       # Need to extract the vector correctly, considering only relevant indices
       final_vector <- as.vector(acc_matrix[,1]) # Full vector matching dataset$mask length

       # Wrap output using create_searchlight_performance (pass the full vector, no IDs needed here)
       out_list[[m_name]] <- create_searchlight_performance(dataset, final_vector, ids=NULL)
   }

   # Create the final searchlight_result object
   structure(
       list(
           results = out_list,
           n_voxels = length(dataset$mask),
           active_voxels = sum(dataset$mask > 0),
           metrics = metric_names
       ),
       class = c("searchlight_result", "list")
   )
}