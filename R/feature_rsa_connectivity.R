#' Extract Per-ROI Predicted and Observed RDM Vectors from Feature RSA Results
#'
#' Convenience helper to pull the compact lower-triangle predicted (and
#' optionally observed) RDM vectors stored by
#' \code{feature_rsa_model(..., return_rdm_vectors = TRUE)} from a
#' \code{regional_mvpa_result}.
#'
#' @param x A \code{regional_mvpa_result} returned by \code{run_regional()} for a
#'   \code{feature_rsa_model}, or a tibble/data frame with columns
#'   \code{roinum} and \code{rdm_vec}.
#'
#' @return A tibble with one row per ROI and columns:
#'   \describe{
#'     \item{roinum}{ROI id.}
#'     \item{n_obs}{Number of observations contributing to the vector.}
#'     \item{observation_index}{List-column of observation ordering used for the
#'       predicted RDM.}
#'     \item{rdm_vec}{List-column containing the lower-triangle predicted RDM
#'       vector for that ROI.}
#'     \item{observed_rdm_vec}{List-column containing the lower-triangle
#'       observed RDM vector for that ROI (if available).}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- run_regional(
#'   feature_rsa_model(dataset, design, method = "pls", return_rdm_vectors = TRUE),
#'   region_mask
#' )
#' vecs <- feature_rsa_rdm_vectors(res)
#' }
#' @export
feature_rsa_rdm_vectors <- function(x) {
  if (is.data.frame(x) && all(c("roinum", "rdm_vec") %in% names(x))) {
    return(tibble::as_tibble(x))
  }

  if (!inherits(x, "regional_mvpa_result")) {
    stop("feature_rsa_rdm_vectors: pass a regional_mvpa_result or a tibble with columns `roinum` and `rdm_vec`.")
  }

  fits <- x$fits
  if (is.null(fits) || !length(fits)) {
    stop("feature_rsa_rdm_vectors: no ROI diagnostics found; re-run feature_rsa_model(..., return_rdm_vectors=TRUE).")
  }

  roi_ids <- seq_along(fits)
  if (!is.null(x$performance_table) &&
      is.data.frame(x$performance_table) &&
      "roinum" %in% names(x$performance_table) &&
      nrow(x$performance_table) == length(fits)) {
    roi_ids <- x$performance_table$roinum
  }

  rows <- lapply(seq_along(fits), function(i) {
    pred <- fits[[i]]
    vec <- tryCatch(pred$predicted_rdm_vec, error = function(...) NULL)
    if (is.null(vec)) {
      return(NULL)
    }
    obs_vec <- tryCatch(pred$observed_rdm_vec, error = function(...) NULL)
    tibble::tibble(
      roinum = roi_ids[[i]],
      n_obs = as.integer(if (is.null(pred$n_obs)) NA_integer_ else pred$n_obs),
      observation_index = list(pred$observation_index),
      rdm_vec = list(as.numeric(vec)),
      observed_rdm_vec = list(if (!is.null(obs_vec)) as.numeric(obs_vec) else NULL)
    )
  })

  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    stop("feature_rsa_rdm_vectors: no predicted RDM vectors found; re-run feature_rsa_model(..., return_rdm_vectors=TRUE).")
  }

  dplyr::bind_rows(rows)
}

.feature_rsa_sparsify_connectivity <- function(mat, keep = 1, absolute = FALSE) {
  if (!is.numeric(keep) || length(keep) != 1L || !is.finite(keep) || keep <= 0 || keep > 1) {
    stop("feature_rsa_connectivity: `keep` must be a scalar in (0, 1].")
  }
  if (keep >= 1 || nrow(mat) < 2L) {
    return(mat)
  }

  out <- mat
  tri <- upper.tri(out)
  vals <- out[tri]
  ok <- is.finite(vals)

  if (!any(ok)) {
    return(out)
  }

  scores <- if (isTRUE(absolute)) abs(vals[ok]) else vals[ok]
  n_keep <- max(1L, ceiling(keep * length(scores)))
  cutoff <- sort(scores, decreasing = TRUE)[n_keep]

  drop_mask <- rep(FALSE, length(vals))
  drop_mask[ok] <- scores < cutoff
  vals[drop_mask] <- 0
  out[tri] <- vals
  out[lower.tri(out)] <- t(out)[lower.tri(out)]
  diag(out) <- 1
  out
}

#' Compute ROI-by-ROI Representational Connectivity from Feature RSA Predictions
#'
#' Forms an ROI x ROI similarity matrix by correlating lower-triangle predicted
#' RDM vectors across ROIs. Sparsification, when requested, is applied only to
#' the final ROI x ROI matrix and never to the per-ROI RDM vectors themselves.
#'
#' @param x Either a \code{regional_mvpa_result} produced by
#'   \code{feature_rsa_model(..., return_rdm_vectors=TRUE)} or the tibble
#'   returned by \code{feature_rsa_rdm_vectors()}.
#' @param method Correlation method used across ROI RDM vectors, one of
#'   \code{"spearman"} or \code{"pearson"}.
#' @param keep Proportion of ROI-ROI edges to retain after optional
#'   sparsification. \code{keep = 1} disables sparsification. For example,
#'   \code{keep = 0.1} retains the top 10\% of finite off-diagonal edges.
#' @param absolute Logical; when \code{TRUE}, rank edges by absolute magnitude
#'   during sparsification. Defaults to \code{FALSE}.
#' @param use Missing-value handling passed to \code{\link[stats]{cor}}.
#'
#' @return A symmetric numeric matrix with ROIs in rows/columns.
#'
#' @examples
#' \dontrun{
#' vecs <- feature_rsa_rdm_vectors(res)
#' conn <- feature_rsa_connectivity(vecs, method = "spearman", keep = 0.1)
#' }
#' @export
feature_rsa_connectivity <- function(x,
                                     method = c("spearman", "pearson"),
                                     keep = 1,
                                     absolute = FALSE,
                                     use = "pairwise.complete.obs") {
  method <- match.arg(method)
  vec_tbl <- feature_rsa_rdm_vectors(x)

  vec_lengths <- vapply(vec_tbl$rdm_vec, length, integer(1))
  if (length(unique(vec_lengths)) != 1L) {
    stop("feature_rsa_connectivity: ROI RDM vectors must all have the same length.")
  }

  obs_idx <- vec_tbl$observation_index
  non_null_obs_idx <- Filter(Negate(is.null), obs_idx)
  if (length(non_null_obs_idx) > 1L) {
    ref_idx <- non_null_obs_idx[[1]]
    same_idx <- vapply(non_null_obs_idx[-1], identical, logical(1), y = ref_idx)
    if (!all(same_idx)) {
      stop("feature_rsa_connectivity: ROI RDM vectors do not share the same observation ordering.")
    }
  }

  vec_mat <- do.call(cbind, vec_tbl$rdm_vec)
  if (!is.matrix(vec_mat)) {
    vec_mat <- matrix(vec_mat, ncol = nrow(vec_tbl))
  }

  conn <- suppressWarnings(stats::cor(vec_mat, method = method, use = use))
  if (!is.matrix(conn)) {
    conn <- matrix(conn, nrow = 1L, ncol = 1L)
  }

  roi_labels <- as.character(vec_tbl$roinum)
  dimnames(conn) <- list(roi_labels, roi_labels)
  diag(conn) <- 1

  .feature_rsa_sparsify_connectivity(conn, keep = keep, absolute = absolute)
}

#' Compute Cross-Connectivity: Predicted-Observed ROI x ROI Matrix
#'
#' Builds an asymmetric ROI x ROI matrix where entry (i, j) is the correlation
#' between the predicted RDM vector of ROI i and the observed RDM vector of
#' ROI j.  This captures how well the model-predicted representational geometry
#' in one ROI matches the data-driven geometry in another.
#'
#' @param x Either a \code{regional_mvpa_result} produced by
#'   \code{feature_rsa_model(..., return_rdm_vectors=TRUE)} or the tibble
#'   returned by \code{feature_rsa_rdm_vectors()}.
#' @param method Correlation method, one of \code{"spearman"} or
#'   \code{"pearson"}.
#' @param use Missing-value handling passed to \code{\link[stats]{cor}}.
#'
#' @return A numeric matrix of dimension n_ROI x n_ROI.  Rows correspond to
#'   predicted RDM vectors and columns to observed RDM vectors.  The matrix is
#'   \emph{not} necessarily symmetric.
#'
#' @examples
#' \dontrun{
#' res <- run_regional(
#'   feature_rsa_model(dataset, design, method = "pls", return_rdm_vectors = TRUE),
#'   region_mask
#' )
#' cross_conn <- feature_rsa_cross_connectivity(res, method = "spearman")
#' }
#' @export
feature_rsa_cross_connectivity <- function(x,
                                           method = c("spearman", "pearson"),
                                           use = "pairwise.complete.obs") {
  method <- match.arg(method)
  vec_tbl <- feature_rsa_rdm_vectors(x)

  ## Check that observed RDM vectors are available
  has_obs <- vapply(vec_tbl$observed_rdm_vec, function(v) !is.null(v), logical(1))
  if (!any(has_obs)) {
    stop("feature_rsa_cross_connectivity: no observed RDM vectors found. ",
         "Re-run with feature_rsa_model(..., return_rdm_vectors = TRUE) ",
         "using a version that stores observed RDM vectors.")
  }
  if (!all(has_obs)) {
    warning("feature_rsa_cross_connectivity: some ROIs lack observed RDM vectors; ",
            "their rows/columns will be NA.")
  }

  ## Validate vector lengths
  pred_lengths <- vapply(vec_tbl$rdm_vec, length, integer(1))
  obs_lengths  <- vapply(vec_tbl$observed_rdm_vec,
                         function(v) if (is.null(v)) NA_integer_ else length(v),
                         integer(1))
  all_lengths <- unique(c(pred_lengths, stats::na.omit(obs_lengths)))
  if (length(all_lengths) != 1L) {
    stop("feature_rsa_cross_connectivity: predicted and observed RDM vectors ",
         "must all have the same length.")
  }

  ## Validate observation ordering
  obs_idx <- vec_tbl$observation_index
  non_null_obs_idx <- Filter(Negate(is.null), obs_idx)
  if (length(non_null_obs_idx) > 1L) {
    ref_idx <- non_null_obs_idx[[1]]
    same_idx <- vapply(non_null_obs_idx[-1], identical, logical(1), y = ref_idx)
    if (!all(same_idx)) {
      stop("feature_rsa_cross_connectivity: ROI RDM vectors do not share ",
           "the same observation ordering.")
    }
  }

  n_roi <- nrow(vec_tbl)
  pred_mat <- do.call(cbind, vec_tbl$rdm_vec)        # n_pairs x n_roi
  obs_mat  <- do.call(cbind, lapply(vec_tbl$observed_rdm_vec, function(v) {
    if (is.null(v)) rep(NA_real_, all_lengths) else v
  }))                                                  # n_pairs x n_roi

  ## Cross-correlation: cor(pred_col_i, obs_col_j)
  cross <- suppressWarnings(stats::cor(pred_mat, obs_mat, method = method, use = use))
  if (!is.matrix(cross)) {
    cross <- matrix(cross, nrow = n_roi, ncol = n_roi)
  }

  roi_labels <- as.character(vec_tbl$roinum)
  dimnames(cross) <- list(predicted = roi_labels, observed = roi_labels)
  cross
}
