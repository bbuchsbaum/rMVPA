# Shared searchlight engine utility functions
# Internal helpers used by both SWIFT and dual LDA fast-path engines.

#' Validate and return combiner function for sampled (randomized/resampled) paths
#' @keywords internal
#' @noRd
.engine_sampled_combiner <- function(combiner, engine_label) {
  if (is.function(combiner)) {
    if (identical(combiner, combine_randomized)) {
      return(combine_randomized)
    }
    stop(
      engine_label, " randomized/resampled fast path currently supports combiner='average' only.",
      call. = FALSE
    )
  }

  if (!is.character(combiner) || length(combiner) == 0L) {
    stop("Invalid randomized/resampled combiner for ", engine_label, " fast path.", call. = FALSE)
  }
  choice <- as.character(combiner)[1]
  if (identical(choice, "average") || identical(choice, "combine_randomized")) {
    return(combine_randomized)
  }
  stop(
    engine_label, " randomized/resampled fast path currently supports combiner='average' only.",
    call. = FALSE
  )
}

#' Extract ROI indices from searchlight iterator
#' @keywords internal
#' @noRd
.engine_extract_roi_indices <- function(slight) {
  lapply(slight, function(roi) {
    ids <- suppressWarnings(as.integer(roi))
    ids <- ids[is.finite(ids)]
    ids <- ids[ids > 0L]
    unique(as.integer(ids))
  })
}

#' Extract ROI center IDs from searchlight iterator
#' @keywords internal
#' @noRd
.engine_extract_roi_ids <- function(slight) {
  if (length(slight) == 0L) {
    return(integer(0))
  }

  if (is.integer(slight[[1]])) {
    ids <- vapply(slight, function(x) {
      val <- attr(x, "center.index")
      if (is.null(val)) NA_integer_ else as.integer(val)[1]
    }, integer(1))
  } else {
    ids <- vapply(slight, function(x) as.integer(x@parent_index), integer(1))
  }
  as.integer(ids)
}

#' Training data of a fast engine: mask columns of `ds$train_data`, checked
#'
#' Rows must match `y` when given. Missing or non-finite values put the run
#' outside the engine's regime.
#' @keywords internal
#' @noRd
.engine_extract_data <- function(ds, engine, mask_indices, y = NULL) {
  x <- as.matrix(neuroim2::series(ds$train_data, mask_indices))
  if (!is.null(y) && nrow(x) != length(y)) {
    stop(engine, ": mismatch between train rows and y_train length.", call. = FALSE)
  }
  if (!all(is.finite(x))) {
    stop(.engine_ineligible(engine, "data contain missing or non-finite values"))
  }
  x
}

#' Mask columns of a sphere's grid coordinates (NA outside the mask)
#' @keywords internal
#' @noRd
.engine_lookup_cols <- function(coords, space_obj, lookup) {
  cc <- lookup[as.integer(neuroim2::grid_to_index(space_obj, coords))]
  cc[cc == 0L] <- NA_integer_
  cc
}

#' Mask columns of every searchlight sphere, streamed once
#'
#' Walks the searchlight once (one sphere at a time, never materialised as a
#' list). With `centers`, only those spheres get columns. Each sphere's
#' columns are in get_searchlight() voxel order, mapped to mask positions
#' through a grid -> column lookup vector built once, so the result is the
#' same as a per-sphere match() against `mask_indices`. Voxels outside the
#' mask come back as NA, as match() would give. With `keep_voxel` (a logical
#' by mask column, the filter_roi() rule), out-of-mask voxels are dropped and
#' a voxel is kept when it passes `keep_voxel` or is the centre itself.
#'
#' @return list(centers = integer centre ids in searchlight order, cols = list
#'   of integer mask-column vectors, one per centre).
#' @keywords internal
#' @noRd
.engine_sphere_cols <- function(ds, radius, centers = NULL, keep_voxel = NULL,
                                mask_indices = NULL, sl = NULL) {
  if (is.null(sl)) sl <- get_searchlight(ds, "standard", radius)
  if (is.null(mask_indices)) {
    mask_indices <- ds$mask_indices %||% compute_mask_indices(ds$mask)
  }
  sp <- neuroim2::space(ds$mask)
  n_space <- prod(dim(sp)[1:3])
  lookup <- integer(n_space)
  lookup[mask_indices] <- seq_along(mask_indices)
  wanted <- NULL
  if (!is.null(centers)) {
    wanted <- logical(n_space)
    wanted[as.integer(centers)] <- TRUE
  }

  n <- length(sl)
  ids <- integer(n)
  cols <- vector("list", n)
  n_kept <- 0L
  for (i in seq_len(n)) {
    w <- sl[[i]]
    pid <- as.integer(w@parent_index)
    if (!is.null(wanted) && !wanted[pid]) next
    cc <- .engine_lookup_cols(w@coords, sp, lookup)
    if (!is.null(keep_voxel)) {
      cc <- cc[!is.na(cc)]
      cc <- cc[keep_voxel[cc] | mask_indices[cc] == pid]
    }
    n_kept <- n_kept + 1L
    ids[n_kept] <- pid
    cols[[n_kept]] <- cc
  }
  list(centers = ids[seq_len(n_kept)], cols = cols[seq_len(n_kept)])
}

#' Mask columns of the given centres, in the order of `centers`
#' @keywords internal
#' @noRd
.engine_sphere_cols_aligned <- function(ds, radius, centers, mask_indices) {
  sc <- .engine_sphere_cols(ds, radius, centers = centers, mask_indices = mask_indices)
  pos <- match(as.integer(centers), sc$centers)
  if (anyNA(pos)) {
    stop("searchlight has no sphere for centre(s): ",
         paste(as.integer(centers)[is.na(pos)], collapse = ", "), call. = FALSE)
  }
  sc$cols[pos]
}

#' filter_roi() rule: a voxel is kept when its range over observations is
#' non-zero (one logical per column of `x`)
#' @keywords internal
#' @noRd
.engine_roi_keep <- function(x) {
  span <- matrixStats::colMaxs(x) - matrixStats::colMins(x)
  is.finite(span) & span > 0
}

#' Grid indices of a failed sphere's voxels, as the general path reports them:
#' voxels passing `keep_voxel`, plus the centre itself
#' @keywords internal
#' @noRd
.engine_roi_indices <- function(mask_indices, cols, centre_id, keep_voxel) {
  cols <- cols[!is.na(cols)]
  as.numeric(mask_indices[cols[keep_voxel[cols] | mask_indices[cols] == centre_id]])
}

#' Failed-centre diagnostics in the layout of the general path's bad_results
#'
#' The general path's error-row columns (see .roi_error_rows()), so fast
#' engines and the general iterator report failed centres identically (the
#' error text, the sphere's indices, and no result).
#' @keywords internal
#' @noRd
.engine_bad_results <- function(ids, error_message, indices = NULL) {
  if (is.null(indices)) indices <- list(integer(0))
  .roi_error_rows(as.integer(ids), error_message, indices = indices)
}

#' Wrap a fast engine's per-centre metrics into a searchlight result
#'
#' `good` marks the centres with a result. Their rows of `perf` go through
#' wrap_out(); `bad_results` (see .engine_bad_results()) is attached to the
#' result, also when no centre is good.
#' @keywords internal
#' @noRd
.engine_wrap <- function(perf, ds, ids, good, bad_results) {
  out <- if (any(good)) {
    wrap_out(perf[good, , drop = FALSE], ds, ids = ids[good])
  } else {
    empty_searchlight_result(ds)
  }
  attr(out, "bad_results") <- bad_results
  out
}
