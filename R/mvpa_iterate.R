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

#' @keywords internal
#' @noRd
.fold_cache_enabled <- function() {
  TRUE
}

#' @keywords internal
#' @noRd
.fold_cache_supported_cv <- function(cv_spec) {
  inherits(cv_spec, c(
    "blocked_cross_validation",
    "kfold_cross_validation",
    "custom_cross_validation"
  ))
}

#' @keywords internal
#' @noRd
.build_cv_fold_cache <- function(mspec) {
  if (is.null(mspec$crossval) || !.fold_cache_supported_cv(mspec$crossval)) {
    return(NULL)
  }
  if (!.fold_cache_enabled()) {
    return(NULL)
  }

  y <- y_train(mspec)
  if (is.null(y)) {
    return(NULL)
  }
  n_obs <- if (is.matrix(y)) nrow(y) else length(y)
  if (is.null(n_obs) || n_obs < 2L) {
    return(NULL)
  }

  index_frame <- tibble::tibble(.row = seq_len(n_obs))
  folds <- crossval_samples(mspec$crossval, index_frame, y)
  if (nrow(folds) == 0L) {
    return(NULL)
  }

  train_idx <- lapply(folds$train, function(x) as.integer(x$idx))
  test_idx <- lapply(folds$test, function(x) as.integer(x$idx))
  fold_id <- if (".id" %in% names(folds)) as.character(folds$.id) else as.character(seq_len(nrow(folds)))

  list(
    cv_class = class(mspec$crossval),
    y = y,
    n_obs = as.integer(n_obs),
    train_idx = train_idx,
    test_idx = test_idx,
    ytrain = lapply(train_idx, function(idx) subset_y(y, idx)),
    ytest = lapply(test_idx, function(idx) subset_y(y, idx)),
    fold_id = fold_id
  )
}

#' @keywords internal
#' @noRd
.resolve_cv_fold_cache <- function(mspec, n_obs) {
  if (!.fold_cache_enabled()) {
    return(NULL)
  }
  cache <- mspec$.cv_fold_cache
  if (is.null(cache) || !is.list(cache)) {
    return(NULL)
  }
  if (is.null(mspec$crossval) || !identical(cache$cv_class, class(mspec$crossval))) {
    return(NULL)
  }
  if (!identical(as.integer(cache$n_obs), as.integer(n_obs))) {
    return(NULL)
  }
  y <- y_train(mspec)
  if (!identical(cache$y, y)) {
    return(NULL)
  }
  if (is.null(cache$train_idx) || is.null(cache$test_idx) || is.null(cache$fold_id)) {
    return(NULL)
  }

  cache
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


#' @keywords internal
#' @noRd
.extract_sample_indices <- function(x) {
  if (is.null(x)) {
    return(integer())
  }
  if (inherits(x, "resample") && !is.null(x$idx)) {
    return(as.integer(x$idx))
  }
  if (is.list(x) && !is.null(x$idx)) {
    return(as.integer(x$idx))
  }
  as.integer(x)
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
external_crossval <- function(mspec, roi, id, center_global_id = NA,
                              xtrain_mat = NULL, xtest_mat = NULL, ind = NULL, ...) {
  matrix_first <- .matrix_first_roi_enabled()

  # Prepare the training data
  xtrain_mat <- if (is.null(xtrain_mat)) {
    as.matrix(neuroim2::values(roi$train_roi))
  } else {
    xtrain_mat
  }
  xtrain <- if (matrix_first) {
    xtrain_mat
  } else {
    tibble::as_tibble(xtrain_mat, .name_repair="minimal")
  }

  ytrain <- y_train(mspec)
 
  # Get the testing labels
  ytest <- y_test(mspec)

  # Get the ROI indices
  if (is.null(ind)) {
    ind <- neuroim2::indices(roi$train_roi)
  }

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
                                          tune_reps=mspec$tune_reps, sl_info = sl_info, quiet_error = TRUE), dots)),
                silent = TRUE)
  
 

  if (inherits(result, "try-error")) {
    # Log a warning if there's an error during model training
    flog.warn("error fitting model %s : %s", id, attr(result, "condition")$message)
    # Store error messages and return a tibble with the error information
    emessage <- if (is.null(attr(result, "condition")$message)) "" else attr(result, "condition")$message
    tibble::tibble(class=list(NULL), probs=list(NULL), y_true=list(ytest),
                   fit=list(NULL), error=TRUE, error_message=emessage)
  } else {
    # Make predictions using the trained model
    xtest_mat <- if (is.null(xtest_mat)) {
      as.matrix(neuroim2::values(roi$test_roi))
    } else {
      xtest_mat
    }
    xtest <- if (matrix_first) {
      xtest_mat
    } else {
      tibble::as_tibble(xtest_mat, .name_repair="minimal")
    }
    pred <- predict(result, xtest, NULL)
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
internal_crossval <- function(mspec, roi, id, center_global_id = NA, x_all = NULL, ind = NULL, ...) {
  matrix_first <- .matrix_first_roi_enabled()
  x_all <- if (is.null(x_all)) as.matrix(neuroim2::values(roi$train_roi)) else x_all
  cv_cache <- .resolve_cv_fold_cache(mspec, n_obs = nrow(x_all))

  # Generate cross-validation samples when no compatible cache is available.
  samples <- if (is.null(cv_cache)) {
    crossval_samples(
      mspec$crossval,
      tibble::as_tibble(x_all, .name_repair = .name_repair),
      y_train(mspec)
    )
  } else {
    tibble::tibble(
      ytrain = cv_cache$ytrain,
      ytest = cv_cache$ytest,
      train = cv_cache$train_idx,
      test = cv_cache$test_idx,
      .id = cv_cache$fold_id
    )
  }

  # Get ROI indices (all global indices in this searchlight/region)
  if (is.null(ind)) {
    ind <- neuroim2::indices(roi$train_roi)
  }
  
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
    train_ind <- .extract_sample_indices(train)
    test_ind <- .extract_sample_indices(test)

    use_matrix_path <- matrix_first && length(train_ind) > 0L && length(test_ind) > 0L
    if (isTRUE(use_matrix_path) || !is.null(cv_cache)) {
      train_mat <- x_all[train_ind, , drop = FALSE]
      test_mat <- x_all[test_ind, , drop = FALSE]
    }

    train_input <- if (isTRUE(use_matrix_path)) {
      train_mat
    } else if (!is.null(cv_cache)) {
      tibble::as_tibble(train_mat, .name_repair=.name_repair)
    } else {
      tibble::as_tibble(train, .name_repair=.name_repair)
    }

    test_input <- if (isTRUE(use_matrix_path)) {
      test_mat
    } else if (!is.null(cv_cache)) {
      tibble::as_tibble(test_mat, .name_repair=.name_repair)
    } else {
      test
    }

    if (length(test_ind) == 0L) {
      test_ind <- as.integer(seq_len(length(ytest)))
    }

    # Check if the number of features is less than 2
    if (ncol(train_input) < 2) {
      # Return an error message
      return(
        format_result(mspec, NULL, error_message="error: less than 2 features",
                      context=list(roi=roi, ytrain=ytrain, ytest=ytest,
                                   train=train_input, test=test_input, test_ind=test_ind, .id=.id))
      )
    }

   
    # Train the model - NOW PASSING sl_info
    tryCatch({
      result <- train_model(mspec,
                            train_input,
                            ytrain,
                            sl_info = sl_info,
                            cv_spec = mspec$crossval,
                            indices = ind,
                            quiet_error = TRUE)
      format_result(mspec, result, error_message = NULL,
                    context = list(roi = roi, ytrain = ytrain, ytest = ytest,
                                   train = train_input, test = test_input,
                                   test_ind = test_ind, .id = .id))
    }, error = function(e) {
      flog.warn("error fitting model %s : %s", id, conditionMessage(e))
      format_result(mspec, result = NULL, error_message = conditionMessage(e),
                    context = list(roi = roi, ytrain = ytrain, ytest = ytest,
                                   train = train_input, test = test_input,
                                   test_ind = test_ind, .id = .id))
    })
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

#' @keywords internal
#' @noRd
.batch_bounds <- function(n_items, batch_size) {
  start <- seq.int(1L, n_items, by = batch_size)
  end <- pmin.int(start + batch_size - 1L, n_items)
  list(start = start, end = end)
}

#' @keywords internal
#' @noRd
.prepare_batch_frame <- function(mod_spec, dset, vox_batch, batch_positions,
                                 batch_rnums, analysis_type, use_shard_backend,
                                 min_voxels_required) {
  size <- as.integer(vapply(vox_batch, length, integer(1)))
  keep_idx <- which(size >= min_voxels_required)
  n_after_size_filter <- length(keep_idx)

  if (n_after_size_filter == 0L) {
    return(list(
      frame = tibble::tibble(),
      n_after_size_filter = 0L,
      get_samples_seconds = 0,
      roi_extract_seconds = 0
    ))
  }

  kept_positions <- as.integer(batch_positions[keep_idx])
  kept_rnums <- batch_rnums[keep_idx]
  kept_sizes <- size[keep_idx]

  t_get_samples <- proc.time()[3]
  sf_samples <- get_samples(mod_spec$dataset, vox_batch[keep_idx])
  get_samples_elapsed <- proc.time()[3] - t_get_samples

  extract_elapsed <- 0
  if (use_shard_backend) {
    sf <- tibble::tibble(
      .id = kept_positions,
      rnum = kept_rnums,
      sample = sf_samples$sample,
      size = kept_sizes
    )
  } else {
    sample_list <- sf_samples$sample
    roi_list <- vector("list", n_after_size_filter)
    t_extract_roi <- proc.time()[3]
    if (analysis_type == "searchlight") {
      for (k in seq_len(n_after_size_filter)) {
        roi_list[[k]] <- extract_roi(
          sample_list[[k]],
          dset,
          center_global_id = kept_rnums[[k]],
          min_voxels = min_voxels_required
        )
      }
    } else {
      for (k in seq_len(n_after_size_filter)) {
        roi_list[[k]] <- extract_roi(
          sample_list[[k]],
          dset,
          min_voxels = min_voxels_required
        )
      }
    }
    extract_elapsed <- proc.time()[3] - t_extract_roi
    sf <- tibble::tibble(
      .id = kept_positions,
      rnum = kept_rnums,
      roi = roi_list,
      size = kept_sizes
    )
  }

  list(
    frame = sf,
    n_after_size_filter = as.integer(n_after_size_filter),
    get_samples_seconds = as.numeric(get_samples_elapsed),
    roi_extract_seconds = as.numeric(extract_elapsed)
  )
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
#' @param drop_probs Logical; if TRUE, drop per-ROI probability matrices after computing metrics. Default FALSE.
#' @param fail_fast Logical; if TRUE, stop immediately on first ROI error. Default FALSE.
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
#' @examples
#' \donttest{
#'   ds <- gen_sample_dataset(c(5,5,5), 20, blocks=2, nlevels=2)
#'   cval <- blocked_cross_validation(ds$design$block_var)
#'   mdl <- load_model("sda_notune")
#'   mspec <- mvpa_model(mdl, ds$dataset, ds$design,
#'     "classification", crossval=cval)
#'   sl <- get_searchlight(ds$dataset, radius=3)
#'   vox_iter <- lapply(sl, function(x) x)
#'   results <- mvpa_iterate(mspec, vox_iter[1:5],
#'     ids=seq_along(vox_iter[1:5]))
#' }
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
  profile_enabled <- isTRUE(getOption("rMVPA.profile_searchlight", FALSE))

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
    batch_bounds <- .batch_bounds(length(ids), batch_size)
    nbatches <- length(batch_bounds$start)
    
    dset <- mod_spec$dataset
    use_shard_backend <- inherits(mod_spec, "shard_model_spec")
    if (!has_test_set(mod_spec)) {
      mod_spec$.cv_fold_cache <- .build_cv_fold_cache(mod_spec)
    }

    if (use_shard_backend) {
      on.exit(shard_cleanup(mod_spec$shard_data), add = TRUE)
    }

    # Precompute nonzero mask indices once per invocation to avoid repeated O(#voxels) scans in as_roi.
    if (!is.null(dset$mask) && is.null(dset$mask_indices)) {
      dset$mask_indices <- compute_mask_indices(dset$mask)
    }
    tot <- length(ids)
    
    results <- vector("list", nbatches)
    skipped_rois <- 0
    processed_rois <- 0
    timing <- if (profile_enabled) {
      list(
        analysis_type = analysis_type,
        backend = if (use_shard_backend) "shard" else "default",
        n_workers = nworkers,
        n_batches = nbatches,
        total_rois = tot,
        batch = vector("list", nbatches),
        totals = list(
          get_samples_seconds = 0,
          roi_extract_seconds = 0,
          run_future_seconds = 0,
          batch_seconds = 0
        )
      )
    } else {
      NULL
    }
    
    for (i in seq_len(nbatches)) {
      tryCatch({
        batch_t0 <- proc.time()[3]
        batch_positions <- seq.int(batch_bounds$start[[i]], batch_bounds$end[[i]])
        batch_rnums <- ids[batch_positions]
        vlist <- vox_list[batch_positions]
        if (verbose) {
          batch_size_current <- length(batch_positions)
          futile.logger::flog.info("Processing batch %s/%s (%s ROIs in this batch)", 
                                  crayon::blue(i), 
                                  crayon::blue(nbatches),
                                  crayon::green(batch_size_current))
        }

        # ---- build sample frame (serial) ----
        futile.logger::flog.debug("Processing batch %d with %d voxels", i, length(vlist))
        # Minimum features allowed (relaxed to 1 for searchlight to keep edge spheres)
        min_voxels_required <- if (analysis_type == "searchlight") 1L else 2L
        batch_prepared <- .prepare_batch_frame(
          mod_spec = mod_spec,
          dset = dset,
          vox_batch = vlist,
          batch_positions = batch_positions,
          batch_rnums = batch_rnums,
          analysis_type = analysis_type,
          use_shard_backend = use_shard_backend,
          min_voxels_required = min_voxels_required
        )
        sf <- batch_prepared$frame
        n_sf_after_size_filter <- batch_prepared$n_after_size_filter
        get_samples_elapsed <- batch_prepared$get_samples_seconds
        extract_elapsed <- batch_prepared$roi_extract_seconds

        futile.logger::flog.debug("Batch %d: get_samples + filter(size>=%d) took %.3f sec",
                                  i, min_voxels_required, get_samples_elapsed)
        
        futile.logger::flog.debug("Sample frame has %d rows after filtering", n_sf_after_size_filter)
        n_sf <- 0L
        run_future_elapsed <- 0
        
        if (n_sf_after_size_filter > 0) {
          if (use_shard_backend) {
            # ---- shard path: skip serial ROI extraction ----
            # Workers will extract ROIs from shared memory in run_future.shard_model_spec.
            # Keep the 'sample' column (data_sample with $vox indices).
            futile.logger::flog.debug("Batch %d: shard backend active, skipping serial ROI extraction", i)
          } else {
            futile.logger::flog.debug("Batch %d: ROI extraction (extract_roi) took %.3f sec",
                                      i, extract_elapsed)
          }

          n_sf <- nrow(sf)
          # ---- parallel processing via run_future ----
          t_run_future <- proc.time()[3]
          results[[i]] <- run_future(mod_spec, sf, processor, verbose,
                                     analysis_type = analysis_type,
                                     drop_probs = drop_probs,
                                     fail_fast = fail_fast)
          run_future_elapsed <- proc.time()[3] - t_run_future
          futile.logger::flog.debug("Batch %d: run_future (parallel section) took %.3f sec",
                                    i, run_future_elapsed)
          processed_rois <- processed_rois + n_sf
          
          futile.logger::flog.debug("Batch %d produced %d results", i, nrow(results[[i]]))

          # Free batch-local ROI data so it doesn't linger until the next
          # iteration reassigns sf.  run_future already freed its copy of
          # the frame; this drops the main-process reference.
          sf <- NULL
          batch_prepared <- NULL
          vlist <- NULL
          gc(FALSE)
        } else {
          skipped_rois <- skipped_rois + length(batch_positions)
          futile.logger::flog.warn("%s Batch %s: All ROIs filtered out (size < 2 voxels)", 
                                  crayon::yellow("!"),
                                  crayon::blue(i))
          results[[i]] <- tibble::tibble(
            result = list(NULL),
            indices = list(NULL),
            performance = list(NULL),
            id = batch_rnums,
            error = TRUE,
            error_message = "ROI filtered out (size < 2 voxels)",
            warning = TRUE,
            warning_message = "ROI filtered out (size < 2 voxels)"
          )
        }
        batch_elapsed <- proc.time()[3] - batch_t0
        if (profile_enabled) {
          timing$totals$get_samples_seconds <- timing$totals$get_samples_seconds + get_samples_elapsed
          timing$totals$roi_extract_seconds <- timing$totals$roi_extract_seconds + extract_elapsed
          timing$totals$run_future_seconds <- timing$totals$run_future_seconds + run_future_elapsed
          timing$totals$batch_seconds <- timing$totals$batch_seconds + batch_elapsed
          timing$batch[[i]] <- list(
            batch_index = i,
            n_rois_requested = length(batch_positions),
            n_rois_after_size_filter = as.integer(n_sf_after_size_filter),
            n_rois_processed = as.integer(n_sf),
            get_samples_seconds = as.numeric(get_samples_elapsed),
            roi_extract_seconds = as.numeric(extract_elapsed),
            run_future_seconds = as.numeric(run_future_elapsed),
            batch_seconds = as.numeric(batch_elapsed)
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
  if (profile_enabled) {
    timing$processed_rois <- processed_rois
    timing$skipped_rois <- skipped_rois
    attr(final_results, "timing") <- timing
  }
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
#' the design's \code{cv_labels} field, set during design construction)
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
  nworkers <- future::nbrOfWorkers()
  chunk_size <- max(1L, ceiling(total_items / (nworkers * 4L)))

  if (is.null(processor)) {
    do_fun <- function(obj, roi, rnum, center_global_id = NA) {
      process_roi(obj, roi, rnum, center_global_id = center_global_id)
    }
    supports_center <- TRUE
  } else {
    do_fun <- processor
    fmls <- tryCatch(names(formals(do_fun)), error = function(...) character(0))
    supports_center <- "center_global_id" %in% fmls || "..." %in% fmls
  }

  if (supports_center) {
    invoke_processor <- function(obj, roi, rnum, center_global_id) {
      do_fun(obj, roi, rnum, center_global_id = center_global_id)
    }
  } else {
    invoke_processor <- function(obj, roi, rnum, center_global_id) {
      do_fun(obj, roi, rnum)
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
        result <- invoke_processor(obj, roi, rnum, center_global_id_to_pass)

        if (!obj$return_predictions) {
          result$result <- list(NULL)
        }
        # Optionally drop dense per-ROI probability matrices after performance
        # has been computed, keeping only prob_observed (one number per trial).
        if (drop_probs) {
          res_list <- result$result
          prob_vals <- vector("list", length(res_list))
          for (j in seq_along(res_list)) {
            res_j <- res_list[[j]]
            if (is.null(res_j) || is.null(res_j$probs)) {
              prob_vals[[j]] <- NULL
              next
            }
            prob_vals[[j]] <- tryCatch(prob_observed(res_j), error = function(e) NULL)
            res_j$probs <- NULL
            res_list[[j]] <- res_j
          }
          result$prob_observed <- prob_vals
          result$result <- res_list
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

  use_progressr <- isTRUE(verbose) &&
    total_items > 0 &&
    requireNamespace("progressr", quietly = TRUE)

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
