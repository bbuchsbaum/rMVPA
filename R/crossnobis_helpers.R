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
#'   All pairwise difference matrices are stacked in a 3-D array and
#'   cross-fold inner products are computed in a vectorised manner using
#'   matrix operations (equivalent to applying `tcrossprod` to each pair's
#'   delta matrix) with diagonals excluded so that only cross-fold terms
#'   contribute.
#'   If `U_folds` contains NAs (e.g., a condition is entirely missing from a fold's
#'   training data), distances involving that condition may result in `NA`.
#'   If you want Mahalanobis distances, ensure `U_folds` contains patterns that have
#'   already been whitened (e.g., by pre-multiplying with \(\Sigma_{noise}^{-1/2}\)).
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

  if (!is.array(U_folds) || length(dim(U_folds)) != 3) {
    rlang::abort("`U_folds` must be a 3D array (Conditions x Voxels x Folds)." )
  }

  K <- dim(U_folds)[1] # Number of conditions
  V <- dim(U_folds)[2] # Number of voxels
  M <- dim(U_folds)[3] # Number of folds

  if (is.null(P_voxels)) {
    P_voxels <- V
  } else {
    if (!is.numeric(P_voxels) || length(P_voxels) != 1 || P_voxels <= 0) {
        rlang::abort("`P_voxels`, if provided, must be a single positive number.")
    }
    if (V != P_voxels) {
      rlang::abort(paste0("Second dimension of `U_folds` (", V,
                         ") must match user-provided `P_voxels` (", P_voxels, ")."))
    }
  }

  if (is.null(dimnames(U_folds)[[1]]) || length(dimnames(U_folds)[[1]]) != K) {
      rlang::abort("First dimension of `U_folds` must have dimnames for conditions.")
  }
  condition_names <- dimnames(U_folds)[[1]]

  if (M < 2) {
    rlang::abort("Crossnobis distance requires at least M=2 folds (partitions). Got M = ", M, ".")
  }
  if (K < 2) {
    if (K == 0) rlang::abort("No conditions found in `U_folds` (dim(U_folds)[1] is 0).")
    # If K=1, no pairs, return empty named vector
    return(setNames(numeric(0), character(0)))
  }

  # Generate all unique unordered pairs of conditions in the exact
  # `lower.tri()` order (row index > column index)
  pair_idx <- which(lower.tri(matrix(TRUE, nrow = K, ncol = K)), arr.ind = TRUE)
  n_pairs <- nrow(pair_idx)
  
  if (n_pairs == 0) { # Should be caught by K < 2, but as a safeguard
      return(setNames(numeric(0), character(0)))
  }

  # Compute delta matrices for all pairs at once (pair x V x M)
  delta_array <- U_folds[pair_idx[, 1], , , drop = FALSE] -
                  U_folds[pair_idx[, 2], , , drop = FALSE]

  # Sum across folds and sum of squares across folds
  delta_sum <- apply(delta_array, c(1, 2), sum)
  delta_sq_sum <- apply(delta_array^2, c(1, 2), sum)

  cross_terms <- delta_sum^2 - delta_sq_sum
  cross_sums <- rowSums(cross_terms)

  crossnobis_distances <- cross_sums / (P_voxels * M * (M - 1))

  pair_names <- paste0(condition_names[pair_idx[, 1]],
                       "_vs_",
                       condition_names[pair_idx[, 2]])

  names(crossnobis_distances) <- pair_names
  return(crossnobis_distances)
} 
