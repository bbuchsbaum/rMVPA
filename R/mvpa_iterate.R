#' @noRd
#' @keywords internal
setup_mvpa_logger <- function() {
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty logging. Please install it.")
  }
  
  # Use the standard layout but with colored messages
  futile.logger::flog.layout(futile.logger::layout.simple)
  
  # Check if we're in a test environment or if user has set a custom threshold
  # Don't override if already set to ERROR or higher (more restrictive)
  current_threshold <- futile.logger::flog.threshold()
  
  # Only set to INFO if current threshold is less restrictive (DEBUG)
  # or if running in interactive mode and not in tests
  if (current_threshold <= futile.logger::DEBUG || 
      (interactive() && Sys.getenv("TESTTHAT") == "")) {
    # Set default threshold to INFO to hide DEBUG messages
    # This prevents common/expected errors from being displayed
    futile.logger::flog.threshold(futile.logger::INFO)
  }
}

#' @keywords internal
try_warning  <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- paste0(warn, str_trim(as.character(w)))
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}



#' @noRd
#' @keywords internal
generate_crossval_samples <- function(mspec, roi) {
  crossval_samples(mspec$crossval, tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair="minimal"), y_train(mspec))
}

#' @noRd
#' @keywords internal
handle_model_training_error <- function(result, id, ytest) {
  futile.logger::flog.warn("Model %s fitting error: %s", 
                          crayon::blue(id), 
                          crayon::red(attr(result, "condition")$message))
  emessage <- if (is.null(attr(result, "condition")$message)) "" else attr(result, "condition")$message
  tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(ytest), 
                 fit=list(NULL), error=TRUE, error_message=emessage)
}

create_result_tibble <- function(cres, ind, mspec, id, result, compute_performance) {
  if (compute_performance) {
    tibble::tibble(result=list(cres), indices=list(ind), 
                   performance=list(compute_performance(mspec, cres)), id=id, 
                   error=FALSE, error_message="~", 
                   warning=!is.null(result$warning), 
                   warning_message=if (is.null(result$warning)) "~" else result$warning)
  } else {
    tibble::tibble(result=list(cres), indices=list(ind), performance=list(NULL), id=id, 
                   error=FALSE, error_message="~", 
                   warning=!is.null(result$warning), 
                   warning_message=if (is.null(result$warning)) "~" else result$warning)
  }
}



#' External Cross-Validation
#'
#' This function performs external cross-validation on the provided ROI and model specification.
#' It returns a tibble with performance metrics, fitted model (optional), and any warnings or errors.
#'
#' @param roi A list containing train_roi and test_roi elements.
#' @param mspec A model specification object.
#' @param id A unique identifier for the model.
#' @param center_global_id Optional global ID of the center voxel. Defaults to NA.
#'
#' @return A tibble with performance metrics, fitted model (optional), and any warnings or errors.
#' @noRd
#' @keywords internal
#' @importFrom stats predict
external_crossval <- function(mspec, roi, id, center_global_id = NA, ...) {
  # Prepare the training data
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair="minimal")

  ytrain <- y_train(mspec)
 
  # Get the testing labels
  ytest <- y_test(mspec)

  # Get the ROI indices
  ind <- neuroim2::indices(roi$train_roi)

  # Determine center_local_id based on center_global_id
  center_local_id <- NA
  if (!is.na(center_global_id)) {
      center_local_id <- match(center_global_id, ind)
      if (is.na(center_local_id)) {
          # Center voxel was filtered out - return graceful skip instead of crashing
          futile.logger::flog.warn("Center voxel %s not present after filtering for ROI %s - skipping", 
                                  center_global_id, id)
          return(tibble::tibble(
            class = list(NULL),
            probs = list(NULL),
            y_true = list(ytest),
            fit = list(NULL),
            error = TRUE,
            error_message = sprintf("Center voxel %s not present after filter_roi; ROI skipped", 
                                  center_global_id),
            warning = TRUE,
            warning_message = "Center voxel removed by filter_roi"
          ))
      }
  }

  # Prepare sl_info
  sl_info <- list(center_local_id = center_local_id, center_global_id = center_global_id)

 
  # Train the model and handle any errors
  # Pass sl_info and other dots
  dots <- list(...)
  result <- try(do.call(train_model, c(list(mspec, xtrain, ytrain, indices=ind, param=mspec$tune_grid, 
                                          tune_reps=mspec$tune_reps, sl_info = sl_info), dots)))
  
 

  if (inherits(result, "try-error")) {
    # Log a warning if there's an error during model training
    flog.warn("error fitting model %s : %s", id, attr(result, "condition")$message)
    # Store error messages and return a tibble with the error information
    emessage <- if (is.null(attr(result, "condition")$message)) "" else attr(result, "condition")$message
    tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(ytest),
                   fit=list(NULL), error=TRUE, error_message=emessage)
  } else {
    # Make predictions using the trained model
    pred <- predict(result, tibble::as_tibble(neuroim2::values(roi$test_roi), .name_repair="minimal"), NULL)
    # Convert predictions to a list
    plist <- lapply(pred, list)
    plist$y_true <- list(ytest)
    plist$test_ind <- list(as.integer(seq_along(ytest)))

    # Create a tibble with the predictions
    ret <- tibble::as_tibble(plist, .name_repair = .name_repair)

    # Wrap the results and return the fitted model if required
    cres <- if (mspec$return_fit) {
      wrap_result(ret, mspec$design, result$fit)
    } else {
      wrap_result(ret, mspec$design)
    }

    # Compute performance and return a tibble with the results and any warnings
    if (mspec$compute_performance) {
      tibble::tibble(result=list(cres), indices=list(ind),
                     performance=list(compute_performance(mspec, cres)), id=id,
                     error=FALSE, error_message="~",
                     warning=!is.null(result$warning),
                     warning_message=if (is.null(result$warning)) "~" else result$warning)
    } else {
      tibble::tibble(result=list(cres), indices=list(ind), performance=list(NULL), id=id,
                     error=FALSE, error_message="~",
                     warning=!is.null(result$warning),
                     warning_message=if (is.null(result$warning)) "~" else result$warning)
    }

  }
}


#' Perform Internal Cross-Validation for MVPA Models
#'
#' This function performs internal cross-validation on a region of interest (ROI) using a specified 
#' MVPA model. It handles the training, prediction, and result aggregation for each cross-validation fold.
#'
#' @param mspec An MVPA model specification object containing:
#'   \describe{
#'     \item{crossval}{Cross-validation specification}
#'     \item{compute_performance}{Logical indicating whether to compute performance metrics}
#'     \item{return_fit}{Logical indicating whether to return fitted models}
#'   }
#' @param roi A list containing at least:
#'   \describe{
#'     \item{train_roi}{Training data as a NeuroVec or NeuroSurfaceVector object}
#'   }
#' @param id Identifier for the current analysis
#' @param center_global_id Optional global ID of the center voxel. Defaults to NA.
#'
#' @return A tibble containing:
#'   \describe{
#'     \item{result}{List of prediction results for each fold}
#'     \item{indices}{ROI indices used in the analysis}
#'     \item{performance}{Performance metrics if compute_performance is TRUE}
#'     \item{id}{Analysis identifier}
#'     \item{error}{Logical indicating if an error occurred}
#'     \item{error_message}{Error message if applicable}
#'     \item{warning}{Logical indicating if a warning occurred}
#'     \item{warning_message}{Warning message if applicable}
#'   }
#'
#' @details
#' The function performs the following steps:
#' 1. Generates cross-validation samples using the specified scheme
#' 2. For each fold:
#'    - Checks for minimum feature requirements
#'    - Trains the model on the training set
#'    - Makes predictions on the test set
#'    - Formats and stores results
#' 3. Merges results across all folds
#'
#' @note
#' This is an internal function used by mvpa_iterate and should not be called directly.
#' It assumes that input validation has already been performed.
#'
#' @keywords internal
#' @noRd
internal_crossval <- function(mspec, roi, id, center_global_id = NA) {
  
  # Generate cross-validation samples
  # Note: This step could potentially be moved outside the function
  samples <- crossval_samples(mspec$crossval, tibble::as_tibble(neuroim2::values(roi$train_roi), 
                                                                .name_repair=.name_repair), y_train(mspec))

  # Get ROI indices (all global indices in this searchlight/region)
  ind <- neuroim2::indices(roi$train_roi)
  
  # Determine the LOCAL index corresponding to the GLOBAL center ID, if provided
  center_local_id <- NA # Default to NA
  if (!is.na(center_global_id)) {
      center_local_id <- match(center_global_id, ind)
      if (is.na(center_local_id)) {
           # Center voxel was filtered out - return graceful skip instead of crashing
           futile.logger::flog.warn("Center voxel %s not present after filtering for ROI %s - skipping", 
                                   center_global_id, id)
           return(tibble::tibble(
             result = list(NULL),
             indices = list(ind),
             performance = list(NULL),
             id = id,
             error = TRUE,
             error_message = sprintf("Center voxel %s not present after filter_roi; ROI skipped", 
                                   center_global_id),
             warning = TRUE,
             warning_message = "Center voxel removed by filter_roi"
           ))
      }
  }

  # Prepare sl_info for train_model (will contain NAs if center_global_id was NA)
  sl_info <- list(
      center_local_id = center_local_id,
      center_global_id = center_global_id 
  )

  # Iterate through the samples and fit the model
  ret <- samples %>% pmap(function(ytrain, ytest, train, test, .id) {
    # Check if the number of features is less than 2
    if (ncol(train) < 2) {
      # Return an error message
      return(
        format_result(mspec, NULL, error_message="error: less than 2 features", context=list(roi=roi, ytrain=ytrain, ytest=ytest, train=train, test=test, .id=.id))
      )
    }

   
    # Train the model - NOW PASSING sl_info
    result <- try(train_model(mspec, 
                              tibble::as_tibble(train, .name_repair=.name_repair), 
                              ytrain,
                              sl_info = sl_info, # Pass sl_info here
                              cv_spec = mspec$crossval, # Pass cv_spec if needed by the train method
                              indices=ind)) # indices might still be useful for some models
   

    # Check if there was an error during model fitting
    if (inherits(result, "try-error")) {
      flog.warn("error fitting model %s : %s", id, attr(result, "condition")$message)
      # Store error messages
      emessage <- if (is.null(attr(result, "condition")$message)) "" else attr(result, "condition")$message
      format_result(mspec, result=NULL, error_message=emessage, context=list(roi=roi, ytrain=ytrain, ytest=ytest, train=train, test=test, .id=.id))
    } else {
      # Predict on test data
      format_result(mspec, result, error_message=NULL, context=list(roi=roi, ytrain=ytrain, ytest=ytest, train=train, test=test, .id=.id))
    }
  }) %>% purrr::discard(is.null) %>% dplyr::bind_rows()
  

  merge_results(mspec, ret, indices=ind, id=id)
}

    


#' @keywords internal
#' @noRd
extract_roi <- function(sample, data, center_global_id = NULL, min_voxels = 2) {
  r <- as_roi(sample,data)

  # Check if as_roi returned an error object (e.g., insufficient voxels after mask filtering)
  if (inherits(r$train_roi, "try-error")) {
    futile.logger::flog.debug("Skipping ROI: as_roi returned error (%s)", as.character(r$train_roi))
    return(NULL)
  }

  # Use silent=TRUE to prevent error messages from being displayed on the console
  # Pass center_global_id to preserve it during filtering (for searchlights)
  # filter_roi will throw an error if < 2 valid voxels remain
  r <- try(filter_roi(r, preserve = center_global_id, min_voxels = min_voxels), silent=TRUE)

  if (inherits(r, "try-error")) {
    futile.logger::flog.debug("Skipping ROI: filter_roi failed (%s)", as.character(r))
    return(NULL)
  }

  r
}
  
#' Iterate MVPA Analysis Over Multiple ROIs
#' 
#' @description
#' Performs multivariate pattern analysis (MVPA) across multiple regions of interest (ROIs) 
#' using batch processing and parallel computation.
#'
#' @param mod_spec An MVPA model specification object containing the dataset to analyze,
#'        compute_performance (logical indicating whether to compute performance metrics),
#'        and return_predictions (logical indicating whether to return predictions).
#' @param vox_list A list of voxel indices or coordinates defining each ROI to analyze.
#' @param ids Vector of identifiers for each ROI analysis. Defaults to 1:length(vox_list).
#' @param batch_size Integer specifying number of ROIs to process per batch.
#'        Defaults to 10\% of total ROIs.
#' @param verbose Logical indicating whether to print progress messages. Defaults to TRUE.
#' @param processor Optional custom processing function. If NULL, uses default processor.
#'        Must accept parameters (obj, roi, rnum) and return a tibble.
#' @param analysis_type Character indicating the type of analysis. Defaults to "searchlight".
#'
#' @details
#' The function processes ROIs in batches to manage memory usage. For each batch:
#' \enumerate{
#'   \item Extracts ROI data from the dataset.
#'   \item Filters out ROIs with fewer than 2 voxels.
#'   \item Processes each ROI using either the default or custom processor.
#'   \item Combines results across all batches.
#' }
#'
#' @return A tibble containing results for each ROI with columns:
#' \describe{
#'   \item{result}{List column of analysis results (NULL if return_predictions=FALSE).}
#'   \item{indices}{List column of ROI indices used.}
#'   \item{performance}{List column of performance metrics (if computed).}
#'   \item{id}{ROI identifier.}
#'   \item{error}{Logical indicating if an error occurred.}
#'   \item{error_message}{Error message if applicable.}
#'   \item{warning}{Logical indicating if a warning occurred.}
#'   \item{warning_message}{Warning message if applicable.}
#' }
#'
#' @importFrom furrr future_pmap
#' @importFrom purrr map pmap
#' @export
mvpa_iterate <- function(mod_spec, vox_list, ids = 1:length(vox_list),
                         batch_size = NULL,
                         verbose = TRUE,
                         processor = NULL,
                         analysis_type = c("searchlight", "regional"),
                         drop_probs = FALSE,
                         fail_fast = FALSE) {
  setup_mvpa_logger()

  if (length(vox_list) == 0) {
    futile.logger::flog.warn("Empty voxel list provided. No analysis to perform.")
    return(tibble::tibble())
  }


  futile.logger::flog.debug("Starting mvpa_iterate with %d voxels", length(vox_list))

  tryCatch({
    assert_that(length(ids) == length(vox_list),
                msg = paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))

    analysis_type <- match.arg(analysis_type)

    # --- auto-size batches based on analysis type and available workers ---
    nworkers <- future::nbrOfWorkers()
    if (is.null(batch_size)) {
      if (analysis_type == "regional") {
        # Regional: few large ROIs. One batch maximises parallelism.
        batch_size <- length(ids)
      } else {
        # Searchlight: many small ROIs. Batch for memory, but keep
        # batches large enough to fill all workers several times over.
        batch_size <- max(nworkers * 20L, as.integer(0.1 * length(ids)))
      }
    }
    batch_size <- max(1L, min(as.integer(batch_size), length(ids)))
    nbatches <- ceiling(length(ids) / batch_size)
    batch_group <- sort(rep(1:nbatches, length.out = length(ids)))
    batch_ids <- split(1:length(ids), batch_group)
    rnums <- split(ids, batch_group)
    
    dset <- mod_spec$dataset
    # Precompute nonzero mask indices once per invocation to avoid repeated O(#voxels) scans in as_roi.
    if (!is.null(dset$mask) && is.null(dset$mask_indices)) {
      dset$mask_indices <- compute_mask_indices(dset$mask)
    }
    tot <- length(ids)
    
    results <- vector("list", length(batch_ids))
    skipped_rois <- 0
    processed_rois <- 0
    
    for (i in seq_along(batch_ids)) {
      tryCatch({
        batch_t0 <- proc.time()[3]
        if (verbose) {
          batch_size_current <- length(batch_ids[[i]])
          futile.logger::flog.info("Processing batch %s/%s (%s ROIs in this batch)", 
                                  crayon::blue(i), 
                                  crayon::blue(nbatches),
                                  crayon::green(batch_size_current))
        }

        # ---- build sample frame (serial) ----
        t_get_samples <- proc.time()[3]
        vlist <- vox_list[batch_ids[[i]]]
        size <- sapply(vlist, function(v) length(v))
        futile.logger::flog.debug("Processing batch %d with %d voxels", i, length(vlist))
        # Minimum features allowed (relaxed to 1 for searchlight to keep edge spheres)
        min_voxels_required <- if (analysis_type == "searchlight") 1L else 2L
        sf <- get_samples(mod_spec$dataset, vox_list[batch_ids[[i]]]) %>% 
          mutate(.id=batch_ids[[i]], rnum=rnums[[i]], size=size) %>% 
          filter(size >= min_voxels_required)
        futile.logger::flog.debug("Batch %d: get_samples + filter(size>=%d) took %.3f sec",
                                  i, min_voxels_required, proc.time()[3] - t_get_samples)
        
        futile.logger::flog.debug("Sample frame has %d rows after filtering", nrow(sf))
        
        if (nrow(sf) > 0) {
          # ---- ROI extraction (serial) ----
          t_extract_roi <- proc.time()[3]
        # For searchlight, pass center_global_id to preserve center during filtering
        if (analysis_type == "searchlight") {
          sf <- sf %>%
            rowwise() %>%
            mutate(roi=list(extract_roi(sample, dset, center_global_id = rnum, min_voxels = min_voxels_required))) %>%
            select(-sample)
        } else {
          sf <- sf %>%
            rowwise() %>%
            mutate(roi=list(extract_roi(sample, dset, min_voxels = min_voxels_required))) %>%
            select(-sample)
        }
          futile.logger::flog.debug("Batch %d: ROI extraction (extract_roi) took %.3f sec",
                                    i, proc.time()[3] - t_extract_roi)

          n_sf <- nrow(sf)
          # ---- parallel processing via run_future ----
          t_run_future <- proc.time()[3]
          results[[i]] <- run_future(mod_spec, sf, processor, verbose,
                                     analysis_type = analysis_type,
                                     drop_probs = drop_probs,
                                     fail_fast = fail_fast)
          futile.logger::flog.debug("Batch %d: run_future (parallel section) took %.3f sec",
                                    i, proc.time()[3] - t_run_future)
          processed_rois <- processed_rois + n_sf
          
          futile.logger::flog.debug("Batch %d produced %d results", i, nrow(results[[i]]))

          # Free batch-local ROI data so it doesn't linger until the next
          # iteration reassigns sf.  run_future already freed its copy of
          # the frame; this drops the main-process reference.
          rm(sf, vlist)
          gc(FALSE)
        } else {
          skipped_rois <- skipped_rois + length(batch_ids[[i]])
          futile.logger::flog.warn("%s Batch %s: All ROIs filtered out (size < 2 voxels)", 
                                  crayon::yellow("!"),
                                  crayon::blue(i))
          results[[i]] <- tibble::tibble(
            result = list(NULL),
            indices = list(NULL),
            performance = list(NULL),
            id = rnums[[i]],
            error = TRUE,
            error_message = "ROI filtered out (size < 2 voxels)",
            warning = TRUE,
            warning_message = "ROI filtered out (size < 2 voxels)"
          )
        }
        
      }, error = function(e) {
        futile.logger::flog.error("Batch %d failed: %s", i, e$message)
        NULL
      })
    }
    
    # Final summary log
    futile.logger::flog.info("\nMVPA Iteration Complete\n- Total ROIs: %s\n- Processed: %s\n- Skipped: %s",
                            crayon::blue(tot),
                            crayon::blue(processed_rois),
                            crayon::yellow(skipped_rois))
  # Combine all results and free the per-batch list immediately
  final_results <- dplyr::bind_rows(results)
  rm(results)
  gc()
  return(final_results)
  }, error = function(e) {
    futile.logger::flog.error("mvpa_iterate failed: %s", e$message)
    # Return a well-formed error tibble so upstream callers don't crash
    tibble::tibble(
      result = list(NULL),
      indices = list(NULL),
      performance = list(NULL),
      id = NA_integer_,
      error = TRUE,
      error_message = paste("mvpa_iterate fatal error:", e$message),
      warning = TRUE,
      warning_message = "mvpa_iterate fatal error"
    )
  })
}

#' Create a lightweight model spec for parallel workers
#'
#' Copies a model specification, excluding the large \code{dataset} field
#' that should never be serialised to parallel workers.
#' ROI data is extracted serially and passed via the sample frame,
#' so workers never need the full dataset.
#'
#' @section Long-term (Tier 3) -- file-backed proxy references:
#' Instead of NULLing the dataset, a future enhancement would replace it
#' with a lightweight proxy holding a file path or \pkg{neuroim2} connection
#' (e.g. \code{H5NeuroVec}).
#' Workers that truly need dataset access (such as
#' \code{y_train.hrfdecoder_model}, which calls \code{nobs(obj$dataset)})
#' could then lazily materialise only the metadata they need, eliminating
#' the serial extraction bottleneck entirely.
#'
#' @param obj A model specification object (inheriting from \code{"model_spec"}).
#' @return A copy of \code{obj} with \code{$dataset} set to \code{NULL},
#'   preserving the S3 class chain for dispatch.
#' @keywords internal
#' @noRd
as_worker_spec <- function(obj) {
  obj$dataset <- NULL
  obj
}

#' @param verbose Logical; print progress messages if \code{TRUE}.
#' @param analysis_type The type of analysis (e.g., "searchlight").
#' @details
#' If the \pkg{progressr} package is installed, this method emits per-task
#' progress updates from parallel workers. Progress handling is enabled for
#' the scope of this call and uses the currently configured handlers.
#' @rdname run_future-methods
#' @export
run_future.default <- function(obj, frame, processor=NULL, verbose=FALSE,
                               analysis_type = "searchlight", drop_probs = FALSE,
                               fail_fast = FALSE, ...) {
  gc()
  # Ensure workers never receive the full dataset.
  obj <- as_worker_spec(obj)
  total_items <- nrow(frame)

  # --- chunk_size controls parallelism and per-worker memory ---
  # Each chunk becomes one future; items within a chunk run serially
  # on a single worker.  We target ~4 chunks per worker for good load
  # balancing while keeping scheduling overhead reasonable.
  # Set options(rMVPA.chunk_size = N) to override:
  #   1   = maximum parallelism, one future per ROI (high overhead)
  #   Inf = original furrr behaviour (fewest futures, highest memory)
  chunk_size <- getOption("rMVPA.chunk_size", NULL)
  if (is.null(chunk_size)) {
    nworkers <- future::nbrOfWorkers()
    chunk_size <- max(1L, ceiling(total_items / (nworkers * 4L)))
  }

  do_fun <- if (is.null(processor)) {
    function(obj, roi, rnum, center_global_id = NA) {
      process_roi(obj, roi, rnum, center_global_id = center_global_id)
    }
  } else {
    processor
  }

  invoke_processor <- function(fun, obj, roi, rnum, center_global_id) {
    fmls <- tryCatch(names(formals(fun)), error = function(...) character(0))
    supports_center <- "center_global_id" %in% fmls || "..." %in% fmls
    if (supports_center) {
      fun(obj, roi, rnum, center_global_id = center_global_id)
    } else {
      fun(obj, roi, rnum)
    }
  }

  run_map <- function(progress_tick = NULL) {
    frame %>% furrr::future_pmap(function(.id, rnum, roi, size) {
      if (!is.null(progress_tick)) {
        on.exit(progress_tick(), add = TRUE)
      }

      tryCatch({
        if (is.null(roi)) {
          # ROI failed validation (e.g. from extract_roi returning NULL due to <2 voxels after filter_roi)
          futile.logger::flog.debug("ROI ID %s: Skipped (failed initial validation in extract_roi, e.g. <2 voxels).", rnum)
          msg <- sprintf("ROI %s failed validation (e.g., <2 voxels after filtering or other extract_roi issue)", rnum)
          if (fail_fast) {
            rlang::abort(message = sprintf("ROI %s failed validation: <2 voxels or invalid after filtering.", rnum))
          }
          return(tibble::tibble(
            result = list(NULL),
            indices = list(NULL),
            performance = list(NULL),
            id = rnum,
            error = TRUE,
            error_message = msg,
            warning = TRUE,
            warning_message = msg
          ))

        }

        # Determine the center_global_id based on analysis type
        center_global_id_to_pass <- if (analysis_type == "searchlight") rnum else NA

        # Pass center_global_id when the processor supports it.
        result <- invoke_processor(do_fun, obj, roi, rnum, center_global_id_to_pass)

        if (!obj$return_predictions) {
          result <- result %>% mutate(result = list(NULL))
        }
        # Optionally drop dense per-ROI probability matrices after performance
        # has been computed, keeping only prob_observed (one number per trial).
        if (drop_probs) {
          result <- result %>%
            dplyr::mutate(
              prob_observed = purrr::map(result, function(res) {
                if (is.null(res)) return(NULL)
                if (!is.null(res$probs)) {
                  tryCatch(prob_observed(res), error = function(e) NULL)
                } else {
                  NULL
                }
              }),
              result = purrr::map(result, function(res) {
                if (is.null(res)) return(NULL)
                if (!is.null(res$probs)) {
                  res$probs <- NULL
                }
                res
              })
            )
        }

        result
      }, error = function(e) {
        # If fail_fast, bubble the error immediately with context + traceback
        if (fail_fast) {
          rlang::abort(
            message = sprintf("ROI %s processing error: %s", rnum, e$message),
            parent = e
          )
        }
        # Capture a short traceback (truncated to avoid memory bloat from
        # thousands of failed edge-of-brain ROIs).
        tb <- tryCatch({
          raw <- paste(utils::capture.output(rlang::trace_back()), collapse = "\n")
          if (nchar(raw) > 500L) paste0(substr(raw, 1L, 500L), "\n... [truncated]") else raw
        }, error = function(...) NA_character_)
        futile.logger::flog.debug("ROI %d: Processing error (%s)", rnum, e$message)
        tibble::tibble(
          result = list(NULL),
          indices = list(NULL),
          performance = list(NULL),
          id = rnum,
          error = TRUE,
          error_message = paste0("Error processing ROI: ", e$message,
                                 if (!is.na(tb)) paste0("\ntrace:\n", tb) else ""),
          warning = TRUE,
          warning_message = paste("Error processing ROI:", e$message)
        )
      })
    }, .options = furrr::furrr_options(seed = TRUE, conditions = "condition",
                                       chunk_size = chunk_size))
  }

  use_progressr <- total_items > 0 && requireNamespace("progressr", quietly = TRUE)

  if (isTRUE(verbose) && total_items > 0 && !use_progressr) {
    futile.logger::flog.debug("progressr unavailable; using batch-level progress logs only.")
  }

  results <- if (use_progressr) {
    progressr::with_progress({
      p <- progressr::progressor(steps = total_items)
      run_map(progress_tick = p)
    }, handlers = progressr::handlers(), enable = TRUE)
  } else {
    run_map(progress_tick = NULL)
  }

  # Free the frame (which holds all batch ROI data) and the closure
  # before GC so the memory is actually reclaimable.
  rm(frame, run_map)
  gc()

  dplyr::bind_rows(results)
}
