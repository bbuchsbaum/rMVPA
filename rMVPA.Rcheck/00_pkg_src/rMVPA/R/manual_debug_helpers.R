# Internal helpers to step through a single searchlight ROI manually

# Safe %||% fallback (used in helpers)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Create an ROI list from a single searchlight window for debugging
#'
#' @keywords internal
as_roi_from_sl <- function(sl_win, dataset, min_voxels = 2, preserve_center = TRUE) {
  # Determine voxel indices from the searchlight window
  vox <- if (!is.null(sl_win@coords)) {
    # neuroim2 ROIVolWindow / ROISurfaceWindow
    sl_win@parent_index
  } else if (!is.null(attr(sl_win, "center.index"))) {
    # Fallback for deflist style with attrs
    attr(sl_win, "center.index")
  } else {
    stop("Unsupported searchlight window object")
  }

  # Mimic data_sample structure used by get_samples / as_roi
  sample_obj <- structure(list(data = NULL, vox = vox), class = "data_sample")

  # Preserve precomputed mask indices if present to avoid O(#voxels) scan
  dset <- dataset
  if (is.null(dset$mask_indices) && !is.null(dset$mask)) {
    dset$mask_indices <- compute_mask_indices(dset$mask)
  }

  extract_roi(sample_obj, dset,
              center_global_id = if (preserve_center) vox else NULL,
              min_voxels = min_voxels)
}

#' Run a single ROI through process_roi for debugging
#'
#' @keywords internal
run_single_roi <- function(mod_spec,
                           sl_win,
                           processor = NULL,
                           preserve_center = TRUE,
                           min_voxels = 2,
                           drop_probs = FALSE) {
  roi <- as_roi_from_sl(sl_win, mod_spec$dataset,
                        min_voxels = min_voxels,
                        preserve_center = preserve_center)

  if (is.null(roi)) {
    return(tibble::tibble(
      result = list(NULL),
      indices = list(NULL),
      performance = list(NULL),
      id = sl_win@parent_index %||% NA_integer_,
      error = TRUE,
      error_message = "ROI failed validation",
      warning = TRUE,
      warning_message = "ROI failed validation",
      context = list(list())
    ))
  }

  proc_fun <- processor %||% function(obj, roi, rnum, ...) {
    process_roi(obj, roi, rnum, ...)
  }

  tryCatch({
    res <- proc_fun(mod_spec, roi, sl_win@parent_index %||% NA_integer_)
    if (!is.data.frame(res)) {
      stop("processor returned non-data.frame result")
    }
    # Optionally drop dense probs
    if (drop_probs && "result" %in% names(res)) {
      res$result <- lapply(res$result, function(r) {
        if (is.null(r)) return(NULL)
        if (!is.null(r$probs)) r$probs <- NULL
        r
      })
    }
    if (!"context" %in% names(res)) {
      res$context <- list(list())
    }
    res
  }, error = function(e) {
    tibble::tibble(
      result = list(NULL),
      indices = list(NULL),
      performance = list(NULL),
      id = sl_win@parent_index %||% NA_integer_,
      error = TRUE,
      error_message = paste("Error processing ROI:", e$message),
      warning = TRUE,
      warning_message = paste("Error processing ROI:", e$message),
      context = list(list())
    )
  })
}
