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
#'   If `U_folds` contains NAs (e.g., a condition is entirely missing from a fold's
#'   training data), distances involving that condition may result in `NA`.
#'   If you want Mahalanobis distances, ensure `U_folds` contains patterns that have
#'   already been whitened (e.g., by pre-multiplying with \(\Sigma_{noise}^{-1/2}\)).
#'
#' @references
#'   Diedrichsen, J., & Kriegeskorte, N. (2017). Representational similarity
#'   analysis – connecting the branches of systems neuroscience.
#'   *Nature Reviews Neuroscience*, 18(3), 179-186.
#'   (Specifically Equation 5 for the unbiased squared Euclidean distance).
#'
#' @importFrom utils combn
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

  # Generate all unique unordered pairs of conditions
  condition_indices_pairs <- utils::combn(K, 2) # 2 x n_pairs matrix
  n_pairs <- ncol(condition_indices_pairs)
  
  if (n_pairs == 0) { # Should be caught by K < 2, but as a safeguard
      return(setNames(numeric(0), character(0)))
  }

  crossnobis_distances <- numeric(n_pairs)
  pair_names <- character(n_pairs)

  for (p_idx in 1:n_pairs) {
    idx_cond1 <- condition_indices_pairs[1, p_idx]
    idx_cond2 <- condition_indices_pairs[2, p_idx]

    # Calculate delta_k_m for all folds m: (U_folds[cond1,,m] - U_folds[cond2,,m])
    # This results in a V x M matrix of difference patterns for the current pair
    # Using drop = FALSE ensures it's always V x M
    delta_patterns_all_folds <- U_folds[idx_cond1, , , drop = FALSE] - U_folds[idx_cond2, , , drop = FALSE]
    
    # Compute matrix of inner products: ip[m,n] = delta_k_m^T %*% delta_k_n
    # t(delta_patterns_all_folds) is M x V
    # delta_patterns_all_folds    is V x M
    # Result is M x M
    inner_product_matrix <- crossprod(delta_patterns_all_folds) # Equivalent to t(X) %*% X for X = delta_patterns_all_folds
    
    # Sum of all m != n terms (diagonal is already zero)
    sum_off_diagonal_ips <- sum(inner_product_matrix, na.rm = FALSE) # Propagate NAs

    # Normalize
    if (M * (M - 1) == 0) { # Should be caught by M < 2
        crossnobis_distances[p_idx] <- NA_real_
    } else {
        crossnobis_distances[p_idx] <- sum_off_diagonal_ips / (P_voxels * M * (M - 1))
    }
    
    pair_names[p_idx] <- paste0(condition_names[idx_cond1], "_vs_", condition_names[idx_cond2])
  }

  names(crossnobis_distances) <- pair_names
  return(crossnobis_distances)
} 