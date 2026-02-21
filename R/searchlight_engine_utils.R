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
