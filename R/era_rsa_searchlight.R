# Dedicated standard-searchlight engine for ERA-RSA.
#
# The general-purpose iterator constructs ROIVec objects for every center, filters those
# objects, converts them back to matrices, and then calls fit_roi(). ERA-RSA
# needs only the filtered train/test matrices and feature indices. This engine
# materializes (or shares) the masked matrices once and preserves the established
# filtering contract before calling the same fit_roi.era_rsa_model() method.

#' @keywords internal
#' @noRd
.is_era_rsa_fast_path <- function(model_spec, method = "standard") {
  dataset <- model_spec$dataset

  inherits(model_spec, "era_rsa_model") &&
    identical(as.character(method)[1L], "standard") &&
    inherits(dataset, "mvpa_image_dataset") &&
    !inherits(dataset, "mvpa_multibasis_image_dataset") &&
    identical(searchlight_scope(dataset), "searchlight") &&
    identical(model_spec$pairing, "one_to_one") &&
    !is.null(dataset$test_data) &&
    isTRUE(has_test_set(model_spec)) &&
    !is.null(output_schema(model_spec))
}

#' @keywords internal
#' @noRd
.resolve_searchlight_engine.era_rsa_model <- function(model_spec, method,
                                                       engine = "auto") {
  requested <- .match_searchlight_engine(engine)
  if (identical(requested, "legacy")) {
    return("legacy")
  }

  eligible <- .is_era_rsa_fast_path(model_spec, method)
  if (identical(requested, "auto")) {
    return(if (eligible) "era_rsa_fast" else "legacy")
  }
  if (identical(requested, "era_rsa_fast") && eligible) {
    return("era_rsa_fast")
  }
  "legacy"
}

#' @keywords internal
#' @noRd
.era_fast_expand_searchlight_ids <- function(roi) {
  ids <- tryCatch({
    coords <- tryCatch(roi@coords, error = function(...) NULL)
    sp <- tryCatch(roi@space, error = function(...) NULL)
    if (!is.null(coords) && !is.null(sp)) {
      neuroim2::grid_to_index(sp, coords)
    } else {
      roi
    }
  }, error = function(...) roi)

  ids <- suppressWarnings(as.integer(ids))
  ids[!is.na(ids) & is.finite(ids) & ids > 0L]
}

#' @keywords internal
#' @noRd
.era_fast_matrix_source <- function(dataset) {
  mask_idx <- dataset$mask_indices
  if (is.null(mask_idx)) {
    mask_idx <- compute_mask_indices(dataset$mask)
  }
  mask_idx <- as.integer(mask_idx)
  if (!length(mask_idx)) {
    stop("ERA-RSA fast searchlight requires at least one active mask voxel.",
         call. = FALSE)
  }

  train_roi <- neuroim2::series_roi(dataset$train_data, mask_idx)
  test_roi <- neuroim2::series_roi(dataset$test_data, mask_idx)
  train <- as.matrix(neuroim2::values(train_roi))
  test <- as.matrix(neuroim2::values(test_roi))

  max_idx <- max(mask_idx)
  spatial_dims <- dim(dataset$mask)
  if (!is.null(spatial_dims) && length(spatial_dims) >= 3L) {
    max_idx <- max(max_idx, prod(spatial_dims[seq_len(3L)]))
  }
  idx_to_col <- integer(max_idx)
  idx_to_col[mask_idx] <- seq_along(mask_idx)

  list(
    type = "matrix",
    train = train,
    test = test,
    idx_to_col = idx_to_col,
    mask_idx = mask_idx
  )
}

#' @keywords internal
#' @noRd
.era_fast_shard_source <- function(shard_data) {
  if (is.null(shard_data$shared_train) || is.null(shard_data$shared_test)) {
    stop("ERA-RSA fast shard searchlight requires shared train and test matrices.",
         call. = FALSE)
  }
  list(
    type = "shard",
    train = shard_data$shared_train,
    test = shard_data$shared_test,
    idx_to_col = shard_data$idx_to_col,
    mask_idx = shard_data$mask_idx
  )
}

#' @keywords internal
#' @noRd
.era_fast_error_result <- function(id, indices = integer(), message) {
  roi_result(
    metrics = NULL,
    indices = indices,
    id = id,
    error = TRUE,
    error_message = message,
    warning = TRUE,
    warning_message = message
  )
}

#' @keywords internal
#' @noRd
.era_fast_filter_features <- function(Xtr, feature_ids, center_id,
                                      min_voxels = 1L) {
  stats <- .filter_roi_stats(Xtr)
  keep <- !stats$nas & stats$sdnonzero

  if (length(center_id) == 1L && !is.na(center_id)) {
    center_col <- match(center_id, feature_ids)
    if (!is.na(center_col) && !keep[center_col]) {
      # This intentionally preserves even an NA or zero-variance center. It is
      # the established searchlight filter_roi() behavior.
      keep[center_col] <- TRUE
    }
  }

  if (sum(keep) < min_voxels) NULL else keep
}

#' @keywords internal
#' @noRd
.era_fast_fit_task <- function(task, source, model, schema_names,
                               fail_fast = FALSE) {
  id <- task$id
  fail <- function(message, indices = integer()) {
    if (isTRUE(fail_fast)) {
      stop(sprintf("ROI %s processing error: %s", id, message), call. = FALSE)
    }
    .era_fast_error_result(id, indices = indices, message = message)
  }

  tryCatch({
    roi_ids <- task$indices
    if (!length(roi_ids)) {
      return(fail("empty ROI"))
    }

    in_bounds <- roi_ids <= length(source$idx_to_col)
    valid <- in_bounds
    valid[in_bounds] <- source$idx_to_col[roi_ids[in_bounds]] > 0L

    # as_roi.data_sample() applies intersect() only when some indices are
    # outside the mask. Reproduce its order and duplicate-removal semantics.
    if (any(!valid)) {
      roi_ids <- unique(roi_ids[valid])
    }
    if (!length(roi_ids)) {
      return(fail("ROI outside mask"))
    }

    cols <- source$idx_to_col[roi_ids]
    Xtr <- source$train[, cols, drop = FALSE]
    keep <- .era_fast_filter_features(
      Xtr = Xtr,
      feature_ids = roi_ids,
      center_id = id,
      min_voxels = 1L
    )
    if (is.null(keep)) {
      return(fail("ROI filtered out", indices = roi_ids))
    }

    kept_ids <- roi_ids[keep]
    roi_data <- list(
      train_data = Xtr[, keep, drop = FALSE],
      test_data = source$test[, cols[keep], drop = FALSE],
      indices = kept_ids,
      train_roi = NULL,
      test_roi = NULL
    )
    context <- list(
      design = model$design,
      cv_spec = model$crossval %||% model$cv_spec,
      id = id,
      center_global_id = id
    )
    res <- fit_roi.era_rsa_model(model, roi_data, context)
    if (!inherits(res, "roi_result")) {
      return(fail("fit_roi.era_rsa_model() did not return an roi_result",
                  indices = kept_ids))
    }
    if (isTRUE(res$error)) {
      if (isTRUE(fail_fast)) {
        stop(res$error_message, call. = FALSE)
      }
      return(res)
    }

    metrics <- unlist(res$metrics, recursive = TRUE, use.names = TRUE)
    if (!is.numeric(metrics) || !identical(names(metrics), schema_names)) {
      return(fail(
        "fit_roi.era_rsa_model() metrics do not match output_schema()",
        indices = kept_ids
      ))
    }
    res$metrics <- metrics
    res
  }, error = function(e) {
    if (isTRUE(fail_fast)) stop(e)
    .era_fast_error_result(
      id,
      message = paste0("Error processing ROI: ", conditionMessage(e))
    )
  })
}

#' @keywords internal
#' @noRd
.era_fast_fit_chunk <- function(tasks, source, model, schema_names,
                                fail_fast = FALSE) {
  lapply(
    tasks,
    .era_fast_fit_task,
    source = source,
    model = model,
    schema_names = schema_names,
    fail_fast = fail_fast
  )
}

#' @keywords internal
#' @noRd
.era_fast_task_chunks <- function(tasks, nchunks) {
  if (!length(tasks)) return(list())
  nchunks <- max(1L, min(as.integer(nchunks), length(tasks)))
  chunk_size <- ceiling(length(tasks) / nchunks)
  bounds <- .batch_bounds(length(tasks), chunk_size)
  Map(
    function(first, last) tasks[first:last],
    bounds$start,
    bounds$end
  )
}

#' @keywords internal
#' @noRd
.era_rsa_standard_searchlight_fast <- function(model_spec, radius, k = NULL,
                                               backend = c("default", "shard", "auto"),
                                               fail_fast = FALSE,
                                               verbose = FALSE) {
  backend <- match.arg(backend)
  if (length(radius) != 1L || is.na(radius) || radius < 1 || radius > 100) {
    stop("ERA-RSA fast searchlight requires one radius between 1 and 100.",
         call. = FALSE)
  }

  runtime_model <- configure_runtime_backend(
    model_spec,
    backend = backend,
    context = "ERA-RSA fast searchlight"
  )
  use_shard_backend <- inherits(runtime_model, "shard_model_spec") &&
    !is.null(runtime_model$shard_data)
  if (use_shard_backend && shard_cleanup_on_exit(runtime_model)) {
    on.exit(shard_cleanup(runtime_model$shard_data), add = TRUE)
  }

  dataset <- runtime_model$dataset
  setup_start <- proc.time()[3]
  slight <- get_searchlight(dataset, "standard", radius, k = k)
  center_ids <- as.integer(get_center_ids(dataset))
  if (length(center_ids) != length(slight)) {
    stop("ERA-RSA fast searchlight: center id count does not match searchlight count.",
         call. = FALSE)
  }

  source <- if (use_shard_backend) {
    .era_fast_shard_source(runtime_model$shard_data)
  } else {
    .era_fast_matrix_source(dataset)
  }
  tasks <- Map(
    function(roi, id) list(
      id = id,
      indices = .era_fast_expand_searchlight_ids(roi)
    ),
    slight,
    center_ids
  )
  rm(slight)

  schema_names <- schema_metric_names(output_schema(runtime_model))
  worker_model <- as_worker_spec(runtime_model)
  worker_model$shard_data <- NULL
  class(worker_model) <- setdiff(class(worker_model), "shard_model_spec")
  setup_elapsed <- proc.time()[3] - setup_start

  run_start <- proc.time()[3]
  nworkers <- future::nbrOfWorkers()
  task_results <- if (use_shard_backend && nworkers > 1L && length(tasks) > 1L) {
    chunks <- .era_fast_task_chunks(tasks, nchunks = nworkers * 4L)
    nested <- future.apply::future_lapply(
      chunks,
      .era_fast_fit_chunk,
      source = source,
      model = worker_model,
      schema_names = schema_names,
      fail_fast = fail_fast,
      future.seed = FALSE,
      future.chunk.size = 1L
    )
    unlist(nested, recursive = FALSE, use.names = FALSE)
  } else {
    .era_fast_fit_chunk(
      tasks,
      source = source,
      model = worker_model,
      schema_names = schema_names,
      fail_fast = fail_fast
    )
  }
  run_elapsed <- proc.time()[3] - run_start

  rows <- dplyr::bind_rows(lapply(task_results, roi_result_to_tibble))
  split <- split_results(rows)
  if (!nrow(split$good)) {
    stop("No valid results for ERA-RSA fast searchlight: all ROIs failed to process.",
         call. = FALSE)
  }
  out <- combine_schema_standard(runtime_model, split$good, split$bad)
  attr(out, "bad_results") <- split$bad
  attr(out, "searchlight_engine") <- "era_rsa_fast"
  if (isTRUE(getOption("rMVPA.profile_searchlight", FALSE))) {
    attr(out, "timing") <- list(
      setup_seconds = as.numeric(setup_elapsed),
      iterate_seconds = as.numeric(run_elapsed),
      combine_seconds = NA_real_,
      n_centers = length(center_ids),
      n_good_results = nrow(split$good),
      n_bad_results = nrow(split$bad),
      backend = if (use_shard_backend) "shard" else "default",
      n_workers = nworkers
    )
  }
  if (isTRUE(verbose)) {
    futile.logger::flog.info(
      "ERA-RSA fast searchlight processed %d centers (%d failed) in %.3f sec",
      length(center_ids), nrow(split$bad), run_elapsed
    )
  }
  out
}

#' @keywords internal
#' @noRd
.run_searchlight_engine.era_rsa_model <- function(model_spec, radius, method,
                                                  engine = "auto",
                                                  niter = 4L,
                                                  combiner = "average",
                                                  drop_probs = FALSE,
                                                  fail_fast = FALSE,
                                                  backend = c("default", "shard", "auto"),
                                                  incremental = TRUE,
                                                  gamma = NULL,
                                                  verbose = FALSE,
                                                  k = NULL,
                                                  ...) {
  requested <- .match_searchlight_engine(engine)
  selected <- .resolve_searchlight_engine.era_rsa_model(
    model_spec,
    method = method,
    engine = requested
  )
  strict_requested <- !identical(requested, "auto") &&
    !identical(requested, "legacy")

  valid_combiner <- identical(combiner, combine_standard) ||
    (is.character(combiner) && length(combiner) > 0L &&
       !is.na(combiner[1L]) && combiner[1L] %in% c("average", "standard"))
  eligible_controls <- !isTRUE(drop_probs) && valid_combiner

  if (!identical(selected, "era_rsa_fast") || !eligible_controls) {
    if (strict_requested) {
      stop(
        sprintf(
          "Requested searchlight engine '%s' is not eligible for this analysis.",
          requested
        ),
        call. = FALSE
      )
    }
    return(list(handled = FALSE, result = NULL, engine = "legacy"))
  }

  fast_res <- try(
    .era_rsa_standard_searchlight_fast(
      model_spec = model_spec,
      radius = radius,
      k = k,
      backend = match.arg(backend),
      fail_fast = fail_fast,
      verbose = verbose
    ),
    silent = TRUE
  )
  if (!inherits(fast_res, "try-error")) {
    return(list(handled = TRUE, result = fast_res, engine = "era_rsa_fast"))
  }

  err_msg <- tryCatch(
    conditionMessage(attr(fast_res, "condition")),
    error = function(...) as.character(fast_res)
  )
  if (strict_requested || isTRUE(fail_fast)) {
    stop(sprintf("Requested searchlight engine 'era_rsa_fast' failed: %s", err_msg),
         call. = FALSE)
  }
  warning(
    sprintf("searchlight engine 'era_rsa_fast' failed, falling back to the general-purpose iterator: %s",
            err_msg),
    call. = FALSE
  )
  list(handled = FALSE, result = NULL, engine = "legacy")
}
