#' Summarize REMAP diagnostics at the ROI level
#'
#' Convenience helper to extract ROI-level REMAP diagnostics from a
#' `regional_mvpa_result`. It prefers metrics recorded in the
#' `performance_table` (so it works even when fits are not returned), and will
#' fall back to averaging diagnostics stored in `fits[[i]]$diag_by_fold` when
#' needed.
#'
#' @param regional_res A `regional_mvpa_result` returned by `run_regional()`.
#'
#' @return A tibble with one row per ROI containing:
#'   - `roinum`: ROI id
#'   - `mean_rank`: mean adapter rank
#'   - `mean_lambda`: mean selected lambda
#'   - `mean_roi_improv`: mean fraction of P→M mismatch explained
#'   - `mean_delta_frob`: mean Frobenius norm of the learned correction
#'
#' @examples
#' \dontrun{
#' res <- run_regional(model, region_mask, return_fits = TRUE)
#' summarize_remap_roi(res)
#' }
#' @export
summarize_remap_roi <- function(regional_res) {
  stopifnot(is.list(regional_res))

  # Prefer performance_table columns if present
  pt <- regional_res$performance_table
  have_pt <- is.data.frame(pt) && all(c("roinum") %in% names(pt))
  if (have_pt) {
    # Gracefully handle missing columns by filling NA
    pick <- function(nm) if (nm %in% names(pt)) pt[[nm]] else rep(NA_real_, nrow(pt))
    out <- tibble::tibble(
      roinum          = pt$roinum,
      mean_rank       = as.numeric(pick("adapter_rank")),
      mean_lambda     = as.numeric(pick("lambda_mean")),
      mean_roi_improv = as.numeric(pick("remap_improv")),
      mean_delta_frob = as.numeric(pick("delta_frob_mean"))
    )
    return(out)
  }

  # Fallback: compute from fits if available
  fits <- regional_res$fits
  if (is.null(fits)) {
    stop("summarize_remap_roi: no performance_table and no fits present; re-run with return_fits=TRUE")
  }

  roi_ids <- vapply(fits, function(x) NA_integer_, FUN.VALUE = integer(1L))
  rows <- lapply(seq_along(fits), function(i) {
    diag <- tryCatch(fits[[i]]$diag_by_fold, error = function(...) NULL)
    if (is.null(diag) || !length(diag)) return(NULL)
    ranks   <- vapply(diag, function(d) d$rank_used,   numeric(1L))
    lambdas <- vapply(diag, function(d) d$lambda_used, numeric(1L))
    improv  <- vapply(diag, function(d) d$roi_improv,  numeric(1L))
    frob    <- vapply(diag, function(d) d$delta_frob,  numeric(1L))
    tibble::tibble(
      roinum          = if (!is.null(regional_res$performance_table$roinum)) regional_res$performance_table$roinum[i] else i,
      mean_rank       = mean(ranks, na.rm = TRUE),
      mean_lambda     = mean(lambdas, na.rm = TRUE),
      mean_roi_improv = mean(pmax(improv, 0), na.rm = TRUE),
      mean_delta_frob = mean(frob, na.rm = TRUE)
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    stop("summarize_remap_roi: no usable diagnostics found.")
  }
  dplyr::bind_rows(rows)
}


#' Summarize REMAP item-level residuals for an ROI
#'
#' Extract per-item distances between memory and perception prototypes in the
#' jointly-whitened space, before (naive) and after the REMAP correction. Works
#' from a `regional_mvpa_result` when `return_fits=TRUE`, or directly from a
#' predictor list returned by `process_roi()`.
#'
#' @param x Either a `regional_mvpa_result` (preferred) or a predictor list
#'   containing `diag_by_fold`.
#' @param roi ROI index or id when `x` is a regional result. If `NULL` and `x`
#'   is a predictor, it is ignored.
#'
#' @return A tibble with columns: `item`, `res_naive`, `res_remap`, `res_ratio`,
#'   and `n_folds`.
#'
#' @examples
#' \dontrun{
#' res <- run_regional(model, region_mask, return_fits = TRUE)
#' items_tbl <- summarize_remap_items(res, roi = 1)
#' }
#' @export
summarize_remap_items <- function(x, roi = NULL) {
  diag_list <- NULL

  if (inherits(x, "regional_mvpa_result")) {
    if (is.null(x$fits)) stop("summarize_remap_items: regional result has no fits; run with return_fits=TRUE or pass a predictor.")
    # resolve roi index
    if (is.null(roi)) stop("summarize_remap_items: please supply roi index/id when passing a regional result.")
    idx <- NULL
    if (length(roi) == 1L && is.numeric(roi)) {
      idx <- as.integer(roi)
    } else {
      # try match by roinum id
      if (!is.null(x$performance_table$roinum)) {
        idx <- match(as.integer(roi), as.integer(x$performance_table$roinum))
      }
    }
    if (is.na(idx) || is.null(idx) || idx < 1L || idx > length(x$fits)) stop("summarize_remap_items: invalid roi index/id.")
    pred <- x$fits[[idx]]
    diag_list <- tryCatch(pred$diag_by_fold, error = function(...) NULL)
  } else if (is.list(x) && !is.null(x$diag_by_fold)) {
    diag_list <- x$diag_by_fold
  } else {
    stop("summarize_remap_items: unsupported input; pass regional result (with return_fits=TRUE) or a predictor list.")
  }

  if (is.null(diag_list) || !length(diag_list)) {
    stop("summarize_remap_items: no diag_by_fold available.")
  }

  rows <- lapply(diag_list, function(d) {
    tibble::tibble(
      item          = as.character(d$train_items),
      res_naive     = as.numeric(d$item_res_naive),
      res_remap     = as.numeric(d$item_res_remap)
    )
  })
  df <- dplyr::bind_rows(rows)
  if (nrow(df) == 0) return(df)
  agg <- df |>
    dplyr::group_by(item) |>
    dplyr::summarise(
      res_naive = mean(res_naive, na.rm = TRUE),
      res_remap = mean(res_remap, na.rm = TRUE),
      n_folds   = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(res_ratio = res_remap / pmax(res_naive, .Machine$double.eps))
  agg
}

