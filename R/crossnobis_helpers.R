#' Compute Crossnobis Gram Components within a Searchlight/ROI
#'
#' Shared implementation for cross-validated second moments and distances.
#'
#' @noRd
.rdm_pair_index_cache <- new.env(parent = emptyenv())

#' @noRd
.rdm_pair_indices <- function(K) {
  if (!is.numeric(K) || length(K) != 1L || is.na(K) || K < 0) {
    rlang::abort("`K` must be a single non-negative integer.")
  }
  K <- as.integer(K)
  key <- as.character(K)
  cached <- .rdm_pair_index_cache[[key]]
  if (!is.null(cached)) {
    return(cached)
  }

  lower_mask <- lower.tri(matrix(TRUE, nrow = K, ncol = K))
  pair_idx <- which(lower_mask, arr.ind = TRUE)
  cached <- list(
    mask = lower_mask,
    pair_idx = pair_idx,
    i = pair_idx[, 1],
    j = pair_idx[, 2]
  )
  .rdm_pair_index_cache[[key]] <- cached
  cached
}

#' @noRd
.crossnobis_pair_indices <- function(K) {
  .rdm_pair_indices(K)
}

#' @noRd
.compute_crossnobis_gram_sl <- function(U_folds) {

  if (!is.array(U_folds) || length(dim(U_folds)) != 3) {
    rlang::abort("`U_folds` must be a 3D array (Conditions x Voxels x Folds)." )
  }
  if (!is.numeric(U_folds)) {
    rlang::abort("`U_folds` must be numeric.")
  }

  K <- dim(U_folds)[1] # Number of conditions
  V <- dim(U_folds)[2] # Number of voxels
  M <- dim(U_folds)[3] # Number of folds

  if (is.null(dimnames(U_folds)[[1]]) || length(dimnames(U_folds)[[1]]) != K) {
      rlang::abort("First dimension of `U_folds` must have dimnames for conditions.")
  }
  condition_names <- dimnames(U_folds)[[1]]

  if (M < 2) {
    rlang::abort("Crossnobis distance requires at least M=2 folds (partitions). Got M = ", M, ".")
  }
  if (K == 0) {
    rlang::abort("No conditions found in `U_folds` (dim(U_folds)[1] is 0).")
  }

  U_sum <- rowSums(U_folds, dims = 2)
  dim(U_sum) <- c(K, V)
  dimnames(U_sum) <- dimnames(U_folds)[1:2]

  sum_gram <- tcrossprod(U_sum)
  U_stacked <- U_folds
  dim(U_stacked) <- c(K, V * M)
  within_fold_gram <- tcrossprod(U_stacked)

  G_cv <- (sum_gram - within_fold_gram) / (M * (M - 1))
  dimnames(G_cv) <- list(condition_names, condition_names)

  list(
    G_cv = G_cv,
    K = K,
    V = V,
    M = M,
    condition_names = condition_names
  )
}

#' Compute Cross-Validated Second Moments within a Searchlight/ROI
#'
#' Calculates the unbiased cross-validated second-moment matrix from independent
#' per-partition condition means.
#'
#' @param U_folds A 3D numeric array (K x V x M) of mean activation patterns.
#' @param vectorize Logical. If TRUE, return the lower-triangle vector in
#'   `lower.tri()` order; otherwise return the full K x K matrix.
#'
#' @return A numeric vector or matrix of cross-validated second moments.
#' @keywords internal
#' @noRd
compute_crossnobis_second_moment_sl <- function(U_folds, vectorize = TRUE) {
  assert_that(is.logical(vectorize), length(vectorize) == 1)
  gram <- .compute_crossnobis_gram_sl(U_folds)
  G_cv <- gram$G_cv

  if (!isTRUE(vectorize)) {
    return(G_cv)
  }

  if (gram$K < 2) {
    return(setNames(numeric(0), character(0)))
  }

  pair_info <- .crossnobis_pair_indices(gram$K)
  out <- G_cv[pair_info$mask]
  names(out) <- paste0(gram$condition_names[pair_info$i],
                       "_vs_",
                       gram$condition_names[pair_info$j])
  out
}

#' Compute Crossnobis Distances within a Searchlight/ROI
#'
#' Calculates the unbiased squared Euclidean distance (Crossnobis distance) for all
#' pairs of conditions based on per-fold mean activation patterns.
#'
#' @param U_folds A 3D numeric array (K x V x M) of mean activation patterns.
#'   K is the number of conditions, V is the number of voxels/features, and
#'   M is the number of cross-validation folds. `dimnames(U_folds)[[1]]`
#'   should provide the condition names.
#' @param P_voxels An integer, the number of voxels/features (V). This is used for
#'   normalization in the Crossnobis formula. If `NULL` (default), it is inferred
#'   from `dim(U_folds)[2]`.
#'
#' @return A named numeric vector of length K*(K-1)/2, containing the Crossnobis
#'   distance for each unique pair of conditions. Names are of the format
#'   "CondA_vs_CondB". The order matches the lower triangle of a distance matrix
#'   vectorized by `lower.tri()`.
#'
#' @details
#' The Crossnobis distance is calculated as:
#'   \[ \tilde d_k = \frac{1}{M(M-1) P} \sum_{m \neq n} \hat\delta_{k,m}^{\top}\, \hat\delta_{k,n} \]
#'   where \(\hat\delta_{k,m} = \hat\mu_{i,m}-\hat\mu_{j,m}\) is the difference pattern for
#'   condition pair \(k=(i,j)\) in fold \(m\), and \(P\) is `P_voxels`.
#'   Cross-fold inner products are computed from Gram matrices, using the
#'   identity \(\sum_{m \neq n} x_m^\top x_n =
#'   \|\sum_m x_m\|^2 - \sum_m \|x_m\|^2\). This avoids materializing a
#'   potentially large condition-pair by voxel by fold array.
#'   If `U_folds` contains NAs (e.g., a condition is entirely missing from a fold's
#'   training data), distances involving that condition may result in `NA`.
#'   If you want Mahalanobis distances, ensure `U_folds` contains patterns that have
#'   already been whitened by W, where W %*% t(W) equals the desired precision matrix.
#'
#' @references
#'   Diedrichsen, J., & Kriegeskorte, N. (2017). Representational similarity
#'   analysis - connecting the branches of systems neuroscience.
#'   *Nature Reviews Neuroscience*, 18(3), 179-186.
#'   (Specifically Equation 5 for the unbiased squared Euclidean distance).
#'
#' @importFrom rlang abort
#' @keywords internal
#' @noRd
compute_crossnobis_distances_sl <- function(U_folds, P_voxels = NULL) {
  gram <- .compute_crossnobis_gram_sl(U_folds)
  K <- gram$K
  V <- gram$V

  if (is.null(P_voxels)) {
    P_voxels <- V
  } else {
    if (!is.numeric(P_voxels) || length(P_voxels) != 1 ||
        !is.finite(P_voxels) || P_voxels <= 0) {
        rlang::abort("`P_voxels`, if provided, must be a single positive number.")
    }
    if (V != P_voxels) {
      rlang::abort(paste0("Second dimension of `U_folds` (", V,
                         ") must match user-provided `P_voxels` (", P_voxels, ")."))
    }
  }

  if (K < 2) {
    # If K=1, no pairs, return empty named vector
    return(setNames(numeric(0), character(0)))
  }

  # Generate all unique unordered pairs of conditions in the exact
  # `lower.tri()` order (row index > column index)
  pair_info <- .crossnobis_pair_indices(K)
  pair_idx <- pair_info$pair_idx
  n_pairs <- nrow(pair_idx)
  
  if (n_pairs == 0) { # Should be caught by K < 2, but as a safeguard
      return(setNames(numeric(0), character(0)))
  }

  G_cv <- gram$G_cv
  pair_i <- pair_info$i
  pair_j <- pair_info$j
  G_diag <- diag(G_cv)
  crossnobis_distances <- (G_diag[pair_i] + G_diag[pair_j] -
    2 * G_cv[cbind(pair_i, pair_j)]) / P_voxels

  pair_names <- paste0(gram$condition_names[pair_i],
                       "_vs_",
                       gram$condition_names[pair_j])

  names(crossnobis_distances) <- pair_names
  return(crossnobis_distances)
}
