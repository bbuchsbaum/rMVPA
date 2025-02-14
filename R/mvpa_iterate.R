#' @noRd
#' @keywords internal
setup_mvpa_logger <- function() {
  if (!requireNamespace("crayon", quietly = TRUE)) {
    stop("Package 'crayon' is required for pretty logging. Please install it.")
  }
  
  # Use the standard layout but with colored messages
  futile.logger::flog.layout(futile.logger::layout.simple)
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
#'
#' @return A tibble with performance metrics, fitted model (optional), and any warnings or errors.
#' @noRd
#' @keywords internal
external_crossval <- function(mspec, roi, id) {
  # Prepare the training data
  xtrain <- tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair="minimal")


  ytrain <- y_train(mspec)
 
  # Get the testing labels
  ytest <- y_test(mspec)

  # Get the ROI indices
  ind <- neuroim2::indices(roi$train_roi)

  #browser()
  # Train the model and handle any errors
  result <- try(train_model(mspec, xtrain, ytrain, indices=ind,
                            param=mspec$tune_grid,
                            tune_reps=mspec$tune_reps))
  
 

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
internal_crossval <- function(mspec, roi, id) {
  # Generate cross-validation samples
  # Note: This step could potentially be moved outside the function
  samples <- crossval_samples(mspec$crossval, tibble::as_tibble(neuroim2::values(roi$train_roi), 
                                                                .name_repair=.name_repair), y_train(mspec))

  # Get ROI indices
  ind <- neuroim2::indices(roi$train_roi)

  # Iterate through the samples and fit the model
  ret <- samples %>% pmap(function(ytrain, ytest, train, test, .id) {
    # Check if the number of features is less than 2
    if (ncol(train) < 2) {
      # Return an error message
      return(
        format_result(mspec, NULL, error_message="error: less than 2 features", context=list(roi=roi, ytrain=ytrain, ytest=ytest, train=train, test=test, .id=.id))
      )
    }

   
   
    # Train the model
    result <- try(train_model(mspec, 
                              tibble::as_tibble(train, .name_repair=.name_repair), 
                              ytrain,
                              indices=ind))
   

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
extract_roi <- function(sample, data) {
  r <- as_roi(sample,data)
  v <- neuroim2::values(r$train_roi)
  r <- try(filter_roi(r))
  if (inherits(r, "try-error") || ncol(v) < 2) {
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
#'        and return_predictions (logical indicating whether to return predictions)
#' @param vox_list A list of voxel indices or coordinates defining each ROI to analyze
#' @param ids Vector of identifiers for each ROI analysis. Defaults to 1:length(vox_list)
#' @param batch_size Integer specifying number of ROIs to process per batch.
#'        Defaults to 10% of total ROIs
#' @param verbose Logical indicating whether to print progress messages. Defaults to TRUE
#' @param processor Optional custom processing function. If NULL, uses default processor.
#'        Must accept parameters (obj, roi, rnum) and return a tibble.
#'
#' @details
#' The function processes ROIs in batches to manage memory usage. For each batch:
#' 1. Extracts ROI data from the dataset
#' 2. Filters out ROIs with fewer than 2 voxels
#' 3. Processes each ROI using either the default or custom processor
#' 4. Combines results across all batches
#'
#' @return A tibble containing results for each ROI with columns:
#' \itemize{
#'   \item{result}{List column of analysis results (NULL if return_predictions=FALSE)}
#'   \item{indices}{List column of ROI indices used}
#'   \item{performance}{List column of performance metrics (if computed)}
#'   \item{id}{ROI identifier}
#'   \item{error}{Logical indicating if an error occurred}
#'   \item{error_message}{Error message if applicable}
#'   \item{warning}{Logical indicating if warning occurred}
#'   \item{warning_message}{Warning message if applicable}
#' }
#'
#' @importFrom furrr future_pmap
#' @importFrom purrr map
#' @export
mvpa_iterate <- function(mod_spec, vox_list, ids=1:length(vox_list), 
                         batch_size=as.integer(.1*length(ids)),
                         verbose=TRUE,
                         processor=NULL) {
  
  setup_mvpa_logger()
  
  if (length(vox_list) == 0) {
    futile.logger::flog.warn("⚠ Empty voxel list provided. No analysis to perform.")
    return(tibble::tibble())
  }
  
  
  
  # Add debugging
  futile.logger::flog.debug("Starting mvpa_iterate with %d voxels", length(vox_list))
  
  tryCatch({
    assert_that(length(ids) == length(vox_list), 
                msg=paste("length(ids) = ", length(ids), "::", "length(vox_list) =", length(vox_list)))
    
    batch_size <- max(1, batch_size)
    nbatches <- ceiling(length(ids)/batch_size)
    batch_group <- sort(rep(1:nbatches, length.out=length(ids)))
    batch_ids <- split(1:length(ids), batch_group)
    rnums <- split(ids, batch_group)
    
    dset <- mod_spec$dataset
    tot <- length(ids)
    
    # Process batches and collect results
    results <- vector("list", length(batch_ids))
    skipped_rois <- 0
    processed_rois <- 0
    
    for(i in seq_along(batch_ids)) {
      tryCatch({
        if(verbose) {
          futile.logger::flog.info("⚡ Processing batch %s/%s", 
                                  crayon::blue(i), 
                                  crayon::blue(nbatches))
        }
        
        vlist <- vox_list[batch_ids[[i]]]
        size <- sapply(vlist, function(v) length(v))
        
        # Add debugging
        futile.logger::flog.debug("Processing batch %d with %d voxels", i, length(vlist))
        #browser()
        sf <- get_samples(mod_spec$dataset, vox_list[batch_ids[[i]]]) %>% 
          mutate(.id=batch_ids[[i]], rnum=rnums[[i]], size=size) %>% 
          filter(size>=2)
        
        # Add debugging
        futile.logger::flog.debug("Sample frame has %d rows after filtering", nrow(sf))
        
        ## gpt01 suggestion::
        #sf <- sf %>%
        #  mutate(roi = map(sample, ~ extract_roi(.x, dset))) %>%
        #  select(-sample)
        
        if (nrow(sf) > 0) {
          sf <- sf %>% 
            rowwise() %>% 
            mutate(roi=list(extract_roi(sample,dset))) %>% 
            select(-sample)
          
          results[[i]] <- run_future(mod_spec, sf, processor, verbose)
          processed_rois <- processed_rois + nrow(sf)
          
          # Add debugging
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
    final_results
  }, error = function(e) {
    futile.logger::flog.error("mvpa_iterate failed: %s", e$message)
    return(tibble::tibble())
  })
}

#' @noRd
run_future.default <- function(obj, frame, processor=NULL, verbose=FALSE, ...) {
  gc()
  total_items <- nrow(frame)
  processed_items <- 0
  
  do_fun <- if (is.null(processor)) {
    function(obj, roi, rnum) {
      process_roi(obj, roi, rnum)
    }
  } else {
    processor
  }
  
  frame %>% furrr::future_pmap(function(.id, rnum, roi, size) {
    # Update progress based on actual items processed
    processed_items <<- processed_items + 1
    if (verbose && (processed_items %% 100 == 0)) {
      progress_percent <- as.integer(processed_items/total_items * 100)
      futile.logger::flog.info("↻ Progress: %s%% complete", 
                              crayon::blue(progress_percent))
    }
    
    result <- do_fun(obj, roi, rnum)
    
    if (!obj$return_predictions) {
      result <- result %>% mutate(result = list(NULL))
    }
    
    result
  }, .options=furrr::furrr_options(seed=TRUE)) %>% 
    purrr::discard(is.null) %>% 
    dplyr::bind_rows()
}








