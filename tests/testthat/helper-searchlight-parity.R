extract_metric_values <- function(x) {
  if (is.null(x)) {
    stop("metric map is NULL", call. = FALSE)
  }
  if (is.atomic(x)) {
    return(as.numeric(x))
  }
  as.numeric(neuroim2::values(x))
}

searchlight_allclose <- function(reference, candidate, atol = 1e-8, rtol = 1e-6) {
  if (length(reference) != length(candidate)) {
    stop("reference/candidate lengths differ", call. = FALSE)
  }
  if (!is.finite(atol) || !is.finite(rtol) || atol < 0 || rtol < 0) {
    stop("atol and rtol must be finite non-negative values", call. = FALSE)
  }

  both_na <- is.na(reference) & is.na(candidate)
  both_pos_inf <- is.infinite(reference) & is.infinite(candidate) & (reference > 0) & (candidate > 0)
  both_neg_inf <- is.infinite(reference) & is.infinite(candidate) & (reference < 0) & (candidate < 0)
  finite_mask <- is.finite(reference) & is.finite(candidate)

  comparable <- both_na | both_pos_inf | both_neg_inf | finite_mask
  non_comparable <- which(!comparable)

  diffs <- abs(reference[finite_mask] - candidate[finite_mask])
  tol <- atol + rtol * pmax(abs(reference[finite_mask]), abs(candidate[finite_mask]))
  bad_finite <- which(diffs > tol)

  list(
    ok = (length(non_comparable) == 0L) && (length(bad_finite) == 0L),
    finite_overlap = sum(finite_mask),
    non_comparable_indices = non_comparable,
    n_bad_finite = length(bad_finite),
    max_abs_diff = if (length(diffs) > 0L) max(diffs) else 0,
    max_tol = if (length(tol) > 0L) max(tol) else 0
  )
}

expect_searchlight_parity <- function(reference,
                                      candidate,
                                      metrics = NULL,
                                      atol = 1e-8,
                                      rtol = 1e-6,
                                      require_overlap = TRUE) {
  ref_metrics <- reference$metrics
  cand_metrics <- candidate$metrics
  if (is.null(ref_metrics) || length(ref_metrics) == 0L) {
    ref_metrics <- names(reference$results)
  }
  if (is.null(cand_metrics) || length(cand_metrics) == 0L) {
    cand_metrics <- names(candidate$results)
  }

  if (is.null(metrics)) {
    expect_setequal(ref_metrics, cand_metrics)
    metrics <- intersect(ref_metrics, cand_metrics)
  } else {
    expect_true(all(metrics %in% ref_metrics))
    expect_true(all(metrics %in% cand_metrics))
  }

  if (!is.null(reference$active_voxels) && !is.null(candidate$active_voxels)) {
    expect_equal(as.integer(reference$active_voxels), as.integer(candidate$active_voxels))
  }

  for (metric in metrics) {
    ref_vals <- extract_metric_values(reference$results[[metric]])
    cand_vals <- extract_metric_values(candidate$results[[metric]])
    expect_equal(length(ref_vals), length(cand_vals))

    cmp <- searchlight_allclose(
      reference = ref_vals,
      candidate = cand_vals,
      atol = atol,
      rtol = rtol
    )
    if (isTRUE(require_overlap)) {
      expect_true(
        cmp$finite_overlap > 0L,
        info = sprintf("metric '%s' has no finite overlap", metric)
      )
    }
    expect_true(
      cmp$ok,
      info = sprintf(
        "metric '%s' parity failed (non-comparable=%d, bad_finite=%d, max_abs_diff=%.6g, max_tol=%.6g)",
        metric,
        length(cmp$non_comparable_indices),
        cmp$n_bad_finite,
        cmp$max_abs_diff,
        cmp$max_tol
      )
    )
  }

  invisible(TRUE)
}
