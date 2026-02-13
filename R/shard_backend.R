#' @title Experimental Shared-Memory Backend for Parallel MVPA
#' @description
#' Uses the \pkg{shard} package to place the full data matrix in POSIX shared
#' memory, enabling zero-copy parallel ROI extraction across \pkg{furrr} workers.
#' Activated per model spec via \code{\link{use_shard}}.
#'
#' @details
#' When enabled, the default serial ROI-extraction loop in
#' \code{\link{mvpa_iterate}} is skipped.
#' Instead, each worker receives lightweight ALTREP handles that point into the
#' shared-memory segment and extracts its own ROI columns on demand.
#' This eliminates the two main memory bottlenecks of the default pipeline:
#' \enumerate{
#'   \item Serial \code{extract_roi()} in the main process.
#'   \item Serialisation of large ROI matrices to workers via \pkg{furrr}.
#' }
#'
#' @section Supported dataset types:
#' \itemize{
#'   \item \code{mvpa_image_dataset} (volumetric)
#'   \item \code{mvpa_surface_dataset} (cortical surface)
#'   \item \code{mvpa_clustered_dataset} (parcellated)
#' }
#'
#' @section Cleanup:
#' Shared-memory segments are released automatically at the end of
#' \code{mvpa_iterate()}.
#' You can also call \code{shard_cleanup(mod_spec$shard_data)} explicitly.
#'
#' @name shard_backend
NULL

# ---- public API ----------------------------------------------------------

#' Enable Shared-Memory Backend for an MVPA Model Spec
#'
#' Prepares the dataset for shared-memory access and tags the model
#' specification so that \code{\link{mvpa_iterate}} uses the shard backend
#' instead of the default \pkg{furrr} pipeline.
#'
#' @param mod_spec A model specification created by \code{\link{mvpa_model}}
#'   (or any constructor that produces an object inheriting from
#'   \code{"model_spec"}).
#'
#' @return A copy of \code{mod_spec} with class \code{"shard_model_spec"}
#'   prepended and a \code{$shard_data} list attached.
#'
#' @details
#' This is an \strong{experimental} feature.
#' It requires the \pkg{shard} package and a platform that supports POSIX
#' shared memory (\code{shm_open}).
#'
#' @examples
#' \dontrun{
#'   mspec <- mvpa_model(mdl, dataset, design, "classification", crossval = cval)
#'   mspec <- use_shard(mspec)
#'   results <- run_searchlight(mspec, ...)
#' }
#'
#' @export
use_shard <- function(mod_spec) {
  if (!requireNamespace("shard", quietly = TRUE)) {
    stop("Package 'shard' is required for the shared-memory backend. ",
         "Install it from source.", call. = FALSE)
  }

  dataset <- mod_spec$dataset
  if (is.null(dataset)) {
    stop("mod_spec$dataset is NULL; cannot prepare shared memory.", call. = FALSE)
  }

  # Clean up old shared handles if use_shard() was called before (audit #4)
  if (!is.null(mod_spec$shard_data)) {
    futile.logger::flog.info("shard backend: cleaning up previous shared-memory handles")
    shard_cleanup(mod_spec$shard_data)
    mod_spec$shard_data <- NULL
  }

  futile.logger::flog.info("shard backend: preparing shared memory for dataset (%s)",
                            paste(class(dataset), collapse = ", "))
  mod_spec$shard_data <- shard_prepare_dataset(dataset)
  class(mod_spec) <- unique(c("shard_model_spec", class(mod_spec)))
  mod_spec
}

#' Resolve Execution Backend for a Model Spec
#'
#' Internal helper used by high-level runners to activate/deactivate the shard
#' backend without requiring explicit \code{use_shard()} calls in user code.
#'
#' @param model_spec A model specification object.
#' @param backend One of \code{"default"}, \code{"shard"}, or \code{"auto"}.
#' @param context Short label used in diagnostic log messages.
#' @return A model specification configured for the selected backend.
#' @keywords internal
#' @noRd
configure_runtime_backend <- function(model_spec,
                                      backend = c("default", "shard", "auto"),
                                      context = "run") {
  backend <- match.arg(backend)

  if (backend == "default") {
    if (inherits(model_spec, "shard_model_spec")) {
      if (!is.null(model_spec$shard_data)) {
        shard_cleanup(model_spec$shard_data)
      }
      model_spec$shard_data <- NULL
      class(model_spec) <- setdiff(class(model_spec), "shard_model_spec")
    }
    return(model_spec)
  }

  if (backend == "shard") {
    if (inherits(model_spec, "shard_model_spec") && !is.null(model_spec$shard_data)) {
      return(model_spec)
    }
    return(use_shard(model_spec))
  }

  # backend == "auto": keep existing shard config if already present
  if (inherits(model_spec, "shard_model_spec") && !is.null(model_spec$shard_data)) {
    return(model_spec)
  }

  # If shard isn't installed, default path
  if (!requireNamespace("shard", quietly = TRUE)) {
    futile.logger::flog.debug(
      "%s backend=auto: package 'shard' unavailable, using default backend", context
    )
    return(model_spec)
  }

  # Try shard and fall back on any setup error (e.g., unsupported dataset type)
  tryCatch(
    use_shard(model_spec),
    error = function(e) {
      futile.logger::flog.debug(
        "%s backend=auto: shard unavailable for this run (%s); using default backend",
        context, e$message
      )
      model_spec
    }
  )
}


# ---- shard_prepare_dataset: S3 generic -----------------------------------

#' Prepare Dataset for Shared-Memory Access
#'
#' S3 generic that extracts the data matrix from a dataset, places it in shared
#' memory via \code{shard::share()}, and builds an index-to-column mapping.
#'
#' @param dataset An \code{mvpa_dataset} (or subclass).
#' @param ... Additional arguments (currently unused).
#' @return A list with shared matrices and metadata including a \code{roi_type}
#'   field (\code{"volumetric"}, \code{"surface"}, or \code{"clustered"}).
#' @keywords internal
#' @noRd
shard_prepare_dataset <- function(dataset, ...) {
  UseMethod("shard_prepare_dataset")
}


#' @keywords internal
#' @noRd
shard_prepare_dataset.mvpa_multibasis_image_dataset <- function(dataset, ...) {
  stop("shard backend: multibasis datasets (mvpa_multibasis_image_dataset) ",
       "are not yet supported. Use the default pipeline for multibasis workflows.",
       call. = FALSE)
}


#' @keywords internal
#' @noRd
shard_prepare_dataset.mvpa_image_dataset <- function(dataset, ...) {
  mask_idx <- dataset$mask_indices
  if (is.null(mask_idx)) {
    mask_idx <- compute_mask_indices(dataset$mask)
  }

  # Extract T x V_masked matrix via series_roi -> values
  full_roi <- neuroim2::series_roi(dataset$train_data, mask_idx)
  train_mat <- neuroim2::values(full_roi)
  rm(full_roi)

  shared_train <- shard::share(train_mat)
  rm(train_mat)

  shared_test <- if (!is.null(dataset$test_data)) {
    test_roi <- neuroim2::series_roi(dataset$test_data, mask_idx)
    test_mat <- neuroim2::values(test_roi)
    rm(test_roi)
    st <- shard::share(test_mat)
    rm(test_mat)
    st
  } else {
    NULL
  }

  # Column mapping: linear voxel index -> column in shared matrix
  max_vox <- prod(dim(dataset$mask)[1:3])
  idx_to_col <- integer(max_vox)
  idx_to_col[mask_idx] <- seq_along(mask_idx)

  brain_space <- neuroim2::space(dataset$mask)

  futile.logger::flog.info(
    "shard backend [volumetric]: shared %d x %d matrix (%d masked voxels)",
    nrow(shared_train), ncol(shared_train), length(mask_idx)
  )

  list(
    roi_type     = "volumetric",
    shared_train = shared_train,
    shared_test  = shared_test,
    idx_to_col   = idx_to_col,
    mask_idx     = mask_idx,
    brain_space  = brain_space
  )
}


#' @keywords internal
#' @noRd
shard_prepare_dataset.mvpa_surface_dataset <- function(dataset, ...) {
  mask_idx <- which(dataset$mask > 0)

  # NeuroSurfaceVector: series_roi returns ROISurfaceVector with @data = T x V
  full_roi <- neuroim2::series_roi(dataset$train_data, mask_idx)
  train_mat <- full_roi@data
  rm(full_roi)

  shared_train <- shard::share(train_mat)
  rm(train_mat)

  shared_test <- if (!is.null(dataset$test_data)) {
    test_roi <- neuroim2::series_roi(dataset$test_data, mask_idx)
    test_mat <- test_roi@data
    rm(test_roi)
    st <- shard::share(test_mat)
    rm(test_mat)
    st
  } else {
    NULL
  }

  # Column mapping: vertex index -> column in shared matrix
  n_nodes <- length(dataset$mask)
  idx_to_col <- integer(n_nodes)
  idx_to_col[mask_idx] <- seq_along(mask_idx)

  geometry <- neurosurf::geometry(dataset$train_data)

  futile.logger::flog.info(
    "shard backend [surface]: shared %d x %d matrix (%d masked vertices)",
    nrow(shared_train), ncol(shared_train), length(mask_idx)
  )

  list(
    roi_type     = "surface",
    shared_train = shared_train,
    shared_test  = shared_test,
    idx_to_col   = idx_to_col,
    mask_idx     = mask_idx,
    geometry     = geometry
  )
}


#' @keywords internal
#' @noRd
shard_prepare_dataset.mvpa_clustered_dataset <- function(dataset, ...) {
  # ClusteredNeuroVec: data lives in @ts (T x K matrix, K = num clusters)
  train_mat <- dataset$train_data@ts

  shared_train <- shard::share(train_mat)
  rm(train_mat)

  shared_test <- if (!is.null(dataset$test_data)) {
    st <- shard::share(dataset$test_data@ts)
    st
  } else {
    NULL
  }

  # For clustered data the "indices" are cluster IDs (1:K).
  # idx_to_col is identity: cluster id -> column
  K <- neuroim2::num_clusters(dataset$cvol)
  idx_to_col <- seq_len(K)

  brain_space <- neuroim2::space(dataset$cvol@mask)
  centroids   <- neuroim2::centroids(dataset$cvol)

  futile.logger::flog.info(
    "shard backend [clustered]: shared %d x %d matrix (%d clusters)",
    nrow(shared_train), ncol(shared_train), K
  )

  list(
    roi_type     = "clustered",
    shared_train = shared_train,
    shared_test  = shared_test,
    idx_to_col   = idx_to_col,
    mask_idx     = seq_len(K),
    brain_space  = brain_space,
    centroids    = centroids
  )
}


#' @keywords internal
#' @noRd
shard_prepare_dataset.default <- function(dataset, ...) {
  stop("shard backend: unsupported dataset class '",
       paste(class(dataset), collapse = "', '"), "'", call. = FALSE)
}


# ---- shard_extract_roi: type-aware ROI construction ----------------------

#' Extract ROI from Shared Memory
#'
#' Called by workers to extract ROI data from the shared-memory matrix and
#' construct proper ROI objects compatible with the standard
#' \code{filter_roi} / \code{process_roi} chain.
#' Dispatches on \code{shard_data$roi_type} to construct the correct
#' ROI class (\code{ROIVec}, \code{ROISurfaceVector}, etc.).
#'
#' @param vox Integer vector of indices for this ROI (linear voxel indices,
#'   vertex indices, or cluster IDs depending on dataset type).
#' @param shard_data List returned by \code{shard_prepare_dataset}.
#' @param center_global_id Optional centre voxel to preserve during filtering.
#' @param min_voxels Minimum features after filtering (default 2).
#' @return A list with \code{train_roi} and \code{test_roi}, or \code{NULL}
#'   if filtering fails.
#' @keywords internal
#' @noRd
shard_extract_roi <- function(vox, shard_data,
                               center_global_id = NULL,
                               min_voxels = 2) {
  switch(shard_data$roi_type,
    volumetric = shard_extract_roi_volumetric(vox, shard_data, center_global_id, min_voxels),
    surface    = shard_extract_roi_surface(vox, shard_data, center_global_id, min_voxels),
    clustered  = shard_extract_roi_clustered(vox, shard_data, center_global_id, min_voxels),
    stop("shard_extract_roi: unknown roi_type '", shard_data$roi_type, "'")
  )
}


#' @keywords internal
#' @noRd
sanitize_shard_indices <- function(vox, max_idx) {
  vox <- as.integer(vox)
  vox[!is.na(vox) & vox > 0L & vox <= max_idx]
}


#' @keywords internal
#' @noRd
shard_extract_roi_volumetric <- function(vox, shard_data,
                                          center_global_id = NULL,
                                          min_voxels = 2) {
  vox <- sanitize_shard_indices(vox, length(shard_data$idx_to_col))
  if (length(vox) < 1L) return(NULL)

  col_idx <- shard_data$idx_to_col[vox]
  valid   <- !is.na(col_idx) & col_idx > 0L
  vox     <- vox[valid]
  col_idx <- col_idx[valid]

  if (length(vox) < 1L) return(NULL)

  train_data <- shard_data$shared_train[, col_idx, drop = FALSE]

  vol_dim <- dim(shard_data$brain_space)
  if (length(vol_dim) > 3L) vol_dim <- vol_dim[1:3]
  coords <- arrayInd(vox, .dim = vol_dim)

  train_roi <- neuroim2::ROIVec(shard_data$brain_space,
                                 coords = coords,
                                 data   = train_data)

  test_roi <- if (!is.null(shard_data$shared_test)) {
    test_data <- shard_data$shared_test[, col_idx, drop = FALSE]
    neuroim2::ROIVec(shard_data$brain_space,
                      coords = coords,
                      data   = test_data)
  } else {
    NULL
  }

  roi <- list(train_roi = train_roi, test_roi = test_roi)
  roi <- try(filter_roi(roi, preserve = center_global_id,
                         min_voxels = min_voxels), silent = TRUE)
  if (inherits(roi, "try-error")) return(NULL)
  roi
}


#' @keywords internal
#' @noRd
shard_extract_roi_surface <- function(vox, shard_data,
                                       center_global_id = NULL,
                                       min_voxels = 2) {
  vox <- sanitize_shard_indices(vox, length(shard_data$idx_to_col))
  if (length(vox) < 1L) return(NULL)

  col_idx <- shard_data$idx_to_col[vox]
  valid   <- !is.na(col_idx) & col_idx > 0L
  vox     <- vox[valid]
  col_idx <- col_idx[valid]

  if (length(vox) < 1L) return(NULL)

  train_data <- shard_data$shared_train[, col_idx, drop = FALSE]

  train_roi <- neurosurf::ROISurfaceVector(
    geometry = shard_data$geometry,
    indices  = vox,
    data     = train_data
  )

  test_roi <- if (!is.null(shard_data$shared_test)) {
    test_data <- shard_data$shared_test[, col_idx, drop = FALSE]
    neurosurf::ROISurfaceVector(
      geometry = shard_data$geometry,
      indices  = vox,
      data     = test_data
    )
  } else {
    NULL
  }

  roi <- list(train_roi = train_roi, test_roi = test_roi)
  roi <- try(filter_roi(roi, preserve = center_global_id,
                         min_voxels = min_voxels), silent = TRUE)
  if (inherits(roi, "try-error")) return(NULL)
  roi
}


#' @keywords internal
#' @noRd
shard_extract_roi_clustered <- function(vox, shard_data,
                                         center_global_id = NULL,
                                         min_voxels = 2) {
  vox <- sanitize_shard_indices(vox, length(shard_data$idx_to_col))
  if (length(vox) < 1L) return(NULL)

  # For clustered data, vox = cluster IDs (neighbors).
  # idx_to_col is identity (1:K -> 1:K columns).
  col_idx <- shard_data$idx_to_col[vox]
  valid   <- !is.na(col_idx) & col_idx > 0L
  vox     <- vox[valid]
  col_idx <- col_idx[valid]

  if (length(vox) < 1L) return(NULL)

  train_data <- shard_data$shared_train[, col_idx, drop = FALSE]

  # Coordinates from cluster centroids (rounded to integer for ROIVec compat)
  roi_coords <- round(shard_data$centroids[vox, , drop = FALSE])

  train_roi <- neuroim2::ROIVec(shard_data$brain_space,
                                 coords = roi_coords,
                                 data   = train_data)

  test_roi <- if (!is.null(shard_data$shared_test)) {
    test_data <- shard_data$shared_test[, col_idx, drop = FALSE]
    neuroim2::ROIVec(shard_data$brain_space,
                      coords = roi_coords,
                      data   = test_data)
  } else {
    NULL
  }

  roi <- list(train_roi = train_roi, test_roi = test_roi)
  roi <- try(filter_roi(roi, preserve = center_global_id,
                         min_voxels = min_voxels), silent = TRUE)
  if (inherits(roi, "try-error")) return(NULL)
  roi
}


# ---- S3 method: run_future.shard_model_spec ------------------------------

#' @rdname run_future-methods
#' @export
run_future.shard_model_spec <- function(obj, frame, processor = NULL,
                                         verbose = FALSE,
                                         analysis_type = "searchlight",
                                         drop_probs = FALSE,
                                         fail_fast = FALSE, ...) {
  gc()

  # Capture shared data before stripping
  shard_data <- obj$shard_data

  # Strip dataset for workers (same as default)
  obj <- as_worker_spec(obj)

  # Remove shard class so downstream process_roi dispatches on the real class
  class(obj) <- setdiff(class(obj), "shard_model_spec")

  total_items <- nrow(frame)

  # Chunk-size heuristic (same as run_future.default)
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

  min_voxels <- if (analysis_type == "searchlight") 1L else 2L

  # ---- parallel map: workers extract ROIs from shared memory -------------
  run_map <- function(progress_tick = NULL) {
    frame %>% furrr::future_pmap(function(.id, rnum, sample, size) {
      if (!is.null(progress_tick)) {
        on.exit(progress_tick(), add = TRUE)
      }

      tryCatch({
        vox <- sample$vox
        center_global_id <- if (analysis_type == "searchlight") rnum else NA

        # ROI extraction from shared memory (zero-copy ALTREP handles)
        roi <- shard_extract_roi(vox, shard_data,
                                  center_global_id = center_global_id,
                                  min_voxels = min_voxels)

        if (is.null(roi)) {
          msg <- sprintf("ROI %s failed validation (shard extract)", rnum)
          if (fail_fast) rlang::abort(msg)
          return(tibble::tibble(
            result = list(NULL), indices = list(NULL),
            performance = list(NULL), id = rnum,
            error = TRUE, error_message = msg,
            warning = TRUE, warning_message = msg
          ))
        }

        center_id <- if (analysis_type == "searchlight") rnum else NA
        result <- invoke_processor(do_fun, obj, roi, rnum, center_id)

        if (!obj$return_predictions) {
          result <- result %>% dplyr::mutate(result = list(NULL))
        }

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
                if (!is.null(res$probs)) res$probs <- NULL
                res
              })
            )
        }

        result
      }, error = function(e) {
        if (fail_fast) {
          rlang::abort(
            message = sprintf("ROI %s processing error: %s", rnum, e$message),
            parent = e
          )
        }
        tb <- tryCatch({
          raw <- paste(utils::capture.output(rlang::trace_back()), collapse = "\n")
          if (nchar(raw) > 500L) {
            paste0(substr(raw, 1L, 500L), "\n... [truncated]")
          } else {
            raw
          }
        }, error = function(...) NA_character_)
        futile.logger::flog.debug("ROI %d: Processing error (%s)", rnum, e$message)
        tibble::tibble(
          result = list(NULL), indices = list(NULL),
          performance = list(NULL), id = rnum,
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

  use_progressr <- total_items > 0 &&
    requireNamespace("progressr", quietly = TRUE)

  if (isTRUE(verbose) && total_items > 0 && !use_progressr) {
    futile.logger::flog.debug(
      "progressr unavailable; using batch-level progress logs only.")
  }

  results <- if (use_progressr) {
    progressr::with_progress({
      p <- progressr::progressor(steps = total_items)
      run_map(progress_tick = p)
    }, handlers = progressr::handlers(), enable = TRUE)
  } else {
    run_map(progress_tick = NULL)
  }

  rm(frame, run_map)
  gc()

  dplyr::bind_rows(results)
}


#' Clean Up Shared-Memory Segments
#'
#' Closes shared-memory segments created by \code{\link{use_shard}}.
#' Called automatically at the end of \code{mvpa_iterate()}; can also be
#' invoked manually.
#'
#' @param shard_data List returned by \code{shard_prepare_dataset}
#'   (i.e.\ \code{mod_spec$shard_data}).
#' @return \code{NULL} (invisibly).
#' @export
shard_cleanup <- function(shard_data) {
  if (is.null(shard_data)) return(invisible(NULL))
  for (nm in c("shared_train", "shared_test")) {
    x <- shard_data[[nm]]
    if (!is.null(x)) {
      tryCatch(close(x), error = function(e) {
        futile.logger::flog.debug(
          "shard_cleanup: close() failed for '%s': %s", nm, e$message)
      })
    }
  }
  invisible(NULL)
}
