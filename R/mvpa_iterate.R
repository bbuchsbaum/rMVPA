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
  futile.logger::flog.warn("⚠ Model %s fitting error: %s", 
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
extract_roi <- function(sample, data, center_global_id = NULL) {
  r <- as_roi(sample,data)
  v <- neuroim2::values(r$train_roi)
  
  # Use silent=TRUE to prevent error messages from being displayed on the console
  # Pass center_global_id to preserve it during filtering (for searchlights)
  r <- try(filter_roi(r, preserve = center_global_id), silent=TRUE)
  
  if (inherits(r, "try-error") || ncol(v) < 2) {
    # Only log at debug level so these expected errors don't alarm users
    if (inherits(r, "try-error")) {
      futile.logger::flog.debug("Skipping ROI: insufficient valid columns")
    } else if (ncol(v) < 2) {
      futile.logger::flog.debug("Skipping ROI: less than 2 columns")
    }
    NULL
  } else {
    r
  }
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
#' @importFrom purrr map
#' @export
mvpa_iterate <- function(mod_spec, vox_list, ids = 1:length(vox_list), 
                         batch_size = as.integer(.1 * length(ids)),
                         verbose = TRUE,
                         processor = NULL,
                         analysis_type = c("searchlight", "regional")) {
  setup_mvpa_logger()
  
  if (length(vox_list) == 0) {
    futile.logger::flog.warn("⚠ Empty voxel list provided. No analysis to perform.")
    return(tibble::tibble())
  }


  
  futile.logger::flog.debug("Starting mvpa_iterate with %d voxels", length(vox_list))
  
  tryCatch({
    assert_that(length(ids) == length(vox_list), 
                msg = paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
    
    analysis_type <- match.arg(analysis_type)

    batch_size <- max(1, batch_size)
    nbatches <- ceiling(length(ids) / batch_size)
    batch_group <- sort(rep(1:nbatches, length.out = length(ids)))
    batch_ids <- split(1:length(ids), batch_group)
    rnums <- split(ids, batch_group)
    
    dset <- mod_spec$dataset
    tot <- length(ids)
    
    results <- vector("list", length(batch_ids))
    skipped_rois <- 0
    processed_rois <- 0
    
    for (i in seq_along(batch_ids)) {
      tryCatch({
        if (verbose) {
          batch_size_current <- length(batch_ids[[i]])
          futile.logger::flog.info("⚡ Processing batch %s/%s (%s ROIs in this batch)", 
                                  crayon::blue(i), 
                                  crayon::blue(nbatches),
                                  crayon::green(batch_size_current))
        }

       
        
        vlist <- vox_list[batch_ids[[i]]]
        size <- sapply(vlist, function(v) length(v))
        
        futile.logger::flog.debug("Processing batch %d with %d voxels", i, length(vlist))
        
        sf <- get_samples(mod_spec$dataset, vox_list[batch_ids[[i]]]) %>% 
          mutate(.id=batch_ids[[i]], rnum=rnums[[i]], size=size) %>% 
          filter(size>=2)
        
        futile.logger::flog.debug("Sample frame has %d rows after filtering", nrow(sf))
        
        if (nrow(sf) > 0) {
          # For searchlight, pass center_global_id to preserve center during filtering
          if (analysis_type == "searchlight") {
            sf <- sf %>% 
              rowwise() %>% 
              mutate(roi=list(extract_roi(sample, dset, center_global_id = rnum))) %>% 
              select(-sample)
          } else {
            sf <- sf %>% 
              rowwise() %>% 
              mutate(roi=list(extract_roi(sample, dset))) %>% 
              select(-sample)
          }
          
          # Strip the dataset from mod_spec before passing to parallel workers
          mod_spec_stripped <- strip_dataset(mod_spec)
          
          # Pass the stripped version and analysis_type
          results[[i]] <- run_future(mod_spec_stripped, sf, processor, verbose,
                                     analysis_type = analysis_type)
          # Clean up any completed futures and run garbage collection between batches
          # Clean up completed futures and run garbage collection
          # Note: ClusterRegistry is no longer exported in future >= 1.34.0
          # Using gc() alone for memory cleanup
          gc()
          processed_rois <- processed_rois + nrow(sf)
          
          futile.logger::flog.debug("Batch %d produced %d results", i, nrow(results[[i]]))
        } else {
          skipped_rois <- skipped_rois + length(batch_ids[[i]])
          futile.logger::flog.warn("%s Batch %s: All ROIs filtered out (size < 2 voxels)", 
                                  crayon::yellow("⚠"),
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
    futile.logger::flog.info("\n✨ MVPA Iteration Complete\n├─ Total ROIs: %s\n├─ Processed: %s\n└─ Skipped: %s",
                            crayon::blue(tot),
                            crayon::blue(processed_rois),
                            crayon::yellow(skipped_rois))
    # Combine all results
    final_results <- dplyr::bind_rows(results)
    # Ensure any remaining futures are cleared and memory is reclaimed
    # Clean up and run garbage collection
    # Note: ClusterRegistry is no longer exported in future >= 1.34.0
    gc()
    return(final_results)
  }, error = function(e) {
    futile.logger::flog.error("mvpa_iterate failed: %s", e$message)
    return(tibble::tibble())
  })
}

#' @param verbose Logical; print progress messages if \code{TRUE}.
#' @param analysis_type The type of analysis (e.g., "searchlight").
#' @rdname run_future-methods
#' @export
run_future.default <- function(obj, frame, processor=NULL, verbose=FALSE, analysis_type, ...) {
  gc()
  total_items <- nrow(frame)
  
  do_fun <- if (is.null(processor)) {
    function(obj, roi, rnum, center_global_id = NA) {
      process_roi(obj, roi, rnum, center_global_id = center_global_id)
    }
  } else {
    processor
  }

  
  results <- frame %>% furrr::future_pmap(function(.id, rnum, roi, size) {
    # Note: Progress tracking removed here because it doesn't work correctly
    # with parallel execution. Each worker would have its own counter,
    # causing duplicate messages and incorrect percentages.
    # Progress is now tracked at the batch level in mvpa_iterate.
    
    tryCatch({
      if (is.null(roi)) {
        # ROI failed validation (e.g. from extract_roi returning NULL due to <2 voxels after filter_roi)
        futile.logger::flog.debug("ROI ID %s: Skipped (failed initial validation in extract_roi, e.g. <2 voxels).", rnum)
        return(tibble::tibble(
          result = list(NULL),
          indices = list(NULL),
          performance = list(NULL),
          id = rnum,
          error = TRUE,
          error_message = "ROI failed validation (e.g., <2 voxels after filtering or other extract_roi issue)",
          warning = TRUE,
          warning_message = "ROI failed validation (e.g., <2 voxels after filtering or other extract_roi issue)"
        ))

      }
      
      # Determine the center_global_id based on analysis type
      center_global_id_to_pass <- if (analysis_type == "searchlight") rnum else NA
      
      # If we are here, ROI is valid. Call the processing function
      # Note: Pass center_global_id_to_pass to the processor if it needs it,
      # otherwise it will be handled if the processor calls internal_crossval
      # The default process_roi likely needs to be adapted or replaced.
      # For now, assuming do_fun (process_roi) will handle it implicitly
      # or pass it down to internal_crossval.
      # TODO: We might need to modify process_roi explicitly or make this call directly
      # to internal_crossval depending on the model type.
      # Let's assume for now the default processor logic implicitly handles internal_crossval
      # and we just need to ensure center_global_id_to_pass is available.
      # If a custom processor is used, it's the user's responsibility to handle it.

      result <- do_fun(obj, roi, rnum, center_global_id = center_global_id_to_pass)
      
      if (!obj$return_predictions) {
        result <- result %>% mutate(result = list(NULL))
      }
      
      result
    }, error = function(e) {
      # Use debug level to avoid alarming users with expected errors
      futile.logger::flog.debug("ROI %d: Processing error (%s)", rnum, e$message)
      tibble::tibble(
        result = list(NULL),
        indices = list(NULL),
        performance = list(NULL),
        id = rnum,
        error = TRUE,
        error_message = paste("Error processing ROI:", e$message),
        warning = TRUE,
        warning_message = paste("Error processing ROI:", e$message)
      )
    })
  }, .options=furrr::furrr_options(seed=TRUE))
  
  # Explicitly cleanup resolved futures to free memory before binding results
  # Final cleanup and garbage collection
  # Note: ClusterRegistry is no longer exported in future >= 1.34.0
  gc()

  results %>% dplyr::bind_rows()
}









