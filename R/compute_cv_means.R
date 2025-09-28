#' Compute Cross-Validated Condition Means within a Searchlight
#'
#' This helper function calculates the mean activation pattern for each condition
#' using data from other cross-validation folds.
#'
#' @param sl_data A numeric matrix (samples x voxels/vertices) containing the data
#'   for the current searchlight.
#' @param mvpa_design The \code{mvpa_design} object associated with the dataset,
#'   containing condition labels and block information.
#' @param cv_spec An object describing the cross-validation scheme, typically created
#'   by functions like \code{\\link{blocked_cross_validation}}, \code{\\link{twofold_blocked_cross_validation}},
#'   \code{\\link{kfold_cross_validation}}, etc. (inheriting from \code{cross_validation}).
#'   This object determines how training/test folds are defined.
#' @param estimation_method Character string specifying the method to estimate means.
#'   Currently supported: \code{"average"} (simple mean of training samples per condition).
#'   \itemize{
#'     \item \code{"average"}: Simple mean of training samples per condition.
#'     \item \code{"L2_norm"}: Identical to \code{"average"} but each condition pattern (row) is finally scaled to unit L2 norm. Useful when you need to equalise overall pattern energy across conditions before RSA.
#'     \item \code{"crossnobis"}: Applies a pre-computed whitening matrix (see `whitening_matrix_W`) to the average pattern of each condition within each cross-validation training fold, before averaging these whitened patterns across folds. This aims to produce noise-normalized condition representations.
#'   }
#'   Default is \code{"average"}.
#' @param whitening_matrix_W Optional V x V numeric matrix, where V is the number of voxels/features in `sl_data`. 
#'   This matrix should be the whitening transformation (e.g., Sigma_noise^(-1/2)) derived from GLM residuals.
#'   Required and used only if `estimation_method = "crossnobis"`.
#' @param return_folds Logical, if TRUE, the function returns a list containing both the
#'   overall mean estimate (`mean_estimate`) and an array of per-fold estimates (`fold_estimates`).
#'   If FALSE (default), only the overall mean estimate is returned.
#'
#' @return If `return_folds = FALSE` (default): A numeric matrix (K x V_sl) where K is the number of 
#'   conditions and V_sl is the number of voxels/vertices in the searchlight. Each row 
#'   represents the cross-validated mean pattern for condition k.
#'   If `return_folds = TRUE`: A list with two elements:
#'   \describe{
#'     \item{`mean_estimate`}{The K x V_sl matrix described above.}
#'     \item{`fold_estimates`}{A K x V_sl x M array, where M is the number of folds,
#'       containing the mean estimate for each condition from each fold.}
#'   }
#'
#' @importFrom stats aggregate
#' @importFrom rlang abort
#' @keywords internal
#' @export
compute_crossvalidated_means_sl <- function(sl_data,
                                              mvpa_design,
                                              cv_spec,
                                              estimation_method = "average",
                                              whitening_matrix_W = NULL,
                                              return_folds = FALSE) {

  # --- Input Checks ---
  if (is.data.frame(sl_data)) {
    sl_data <- as.matrix(sl_data)
  } else if (inherits(sl_data, "Matrix")) {
    sl_data <- as.matrix(sl_data)
  }
  if (is.matrix(sl_data) && !is.numeric(sl_data)) {
    storage.mode(sl_data) <- "double"
  }
  if (!is.matrix(sl_data) || !is.numeric(sl_data)) {
    rlang::abort("`sl_data` must be a numeric matrix (samples x features).")
  }
  if (!inherits(mvpa_design, "mvpa_design")) {
    rlang::abort("`mvpa_design` must be an object of class 'mvpa_design'.")
  }
  
  if (!inherits(cv_spec, "cross_validation")) {
    rlang::abort("`cv_spec` must be an object inheriting from class 'cross_validation'.")
  }
  design_rows <- NULL
  if (!is.null(mvpa_design$design_matrix)) {
    design_rows <- nrow(mvpa_design$design_matrix)
  } else if (!is.null(mvpa_design$train_design)) {
    design_rows <- nrow(mvpa_design$train_design)
  } else if (!is.null(mvpa_design$Y)) {
    design_rows <- length(mvpa_design$Y)
  } else if (!is.null(mvpa_design$y_train)) {
    design_rows <- length(mvpa_design$y_train)
  }
  if (!is.null(design_rows) && nrow(sl_data) != design_rows) {
     # Assuming design-related rows correspond to samples
     rlang::abort("Number of rows in `sl_data` must match samples implied by `mvpa_design`.")
  }

  estimation_method <- match.arg(estimation_method, choices = c("average", "L2_norm", "crossnobis"))

  if (estimation_method == "crossnobis") {
    if (is.null(whitening_matrix_W)) {
      rlang::abort("`whitening_matrix_W` must be provided when `estimation_method = crossnobis`.")
    }
    if (!is.matrix(whitening_matrix_W) || !is.numeric(whitening_matrix_W)) {
      rlang::abort("`whitening_matrix_W` must be a numeric matrix.")
    }
    V_sl <- ncol(sl_data)
    if (nrow(whitening_matrix_W) != V_sl || ncol(whitening_matrix_W) != V_sl) {
      rlang::abort(paste0("`whitening_matrix_W` must be a square matrix with dimensions V_sl x V_sl (",
                        V_sl, "x", V_sl, "), matching the number of features in `sl_data`."))
    }
  }

  assert_that(is.logical(return_folds), length(return_folds) == 1)

  # --- Preparation ---
  condition_labels <- mvpa_design$Y # Assuming Y holds condition labels
  if (is.null(condition_labels)) {
      condition_labels <- mvpa_design$y_train
  }
  if (is.null(condition_labels)) {
      condition_labels <- y_train(mvpa_design)
  }
  if (is.null(condition_labels)) {
      rlang::abort("`mvpa_design` must contain condition labels, typically in `$Y`.")
  }
  # Preserve a consistent condition order. Prefer factor levels if `condition_labels` is a factor,
  # otherwise use the order of appearance.
  if (is.factor(condition_labels) && !is.null(levels(condition_labels))) {
      unique_conditions <- levels(condition_labels)
  } else {
      unique_conditions <- unique(condition_labels)
  }
  n_conditions <- length(unique_conditions)
  n_voxels <- ncol(sl_data)
  n_folds <- get_nfolds(cv_spec)

  # Confirm the train_indices generic is available before using it
  if (!exists("train_indices", mode = "function")) {
    rlang::abort("Required generic 'train_indices' is not available in scope.")
  }

  # Pre-compute training indices for each fold and label factor once
  train_idx_list <- lapply(seq_len(n_folds), function(f) train_indices(cv_spec, f))
  condition_labels_factor <- factor(condition_labels, levels = unique_conditions)

  # Prepare containers for online averaging (memory-efficient vs full 3-D array)
  U_hat_sl_cum <- matrix(0, nrow = n_conditions, ncol = n_voxels,
                         dimnames = list(unique_conditions, colnames(sl_data)))
  cond_fold_counts <- integer(n_conditions) # how many folds contributed to each condition

  # Initialize array for per-fold estimates if requested
  if (return_folds) {
    U_folds_array <- array(NA_real_,
                           dim = c(n_conditions, n_voxels, n_folds),
                           dimnames = list(unique_conditions, colnames(sl_data), paste0("Fold", 1:n_folds)))
  } else {
    U_folds_array <- NULL # Not used, keep it NULL
  }

  # --- Calculate Means per Fold ---
  # Iterate folds and accumulate means using vectorised operations
  for (i in seq_len(n_folds)) {

    fold_train_idx <- train_idx_list[[i]]

    if (length(fold_train_idx) == 0) {
      warning(paste("Fold", i, "has no training samples. Skipping."))
      if (return_folds) {
        U_folds_array[, , i] <- NA_real_
      }
      next
    }

    train_data_fold <- sl_data[fold_train_idx, , drop = FALSE]
    train_labels_fold <- condition_labels_factor[fold_train_idx]

    # Compute sums and counts per condition efficiently
    sums_by_cond <- rowsum(train_data_fold, group = train_labels_fold,
                           reorder = FALSE, na.rm = TRUE)
    counts_by_cond <- tabulate(train_labels_fold, nbins = n_conditions)

    fold_means <- matrix(NA_real_, nrow = n_conditions, ncol = n_voxels,
                         dimnames = list(unique_conditions, colnames(sl_data)))
    present_idx <- which(counts_by_cond > 0)
    if (length(present_idx) > 0) {
      aligned_sums <- matrix(0, nrow = n_conditions, ncol = n_voxels,
                             dimnames = list(unique_conditions, colnames(sl_data)))
      if (nrow(sums_by_cond) > 0) {
        sums_by_cond_mat <- as.matrix(sums_by_cond)
        match_idx <- match(rownames(sums_by_cond_mat), unique_conditions)
        valid_match <- !is.na(match_idx)
        if (any(!valid_match)) {
          warning("Encountered training condition without a matching level in unique_conditions when aggregating fold means. Ignoring unmatched rows.")
        }
        if (any(valid_match)) {
          aligned_sums[match_idx[valid_match], ] <- sums_by_cond_mat[valid_match, , drop = FALSE]
        }
      }
      fold_means[present_idx, ] <- aligned_sums[present_idx, , drop = FALSE] /
        counts_by_cond[present_idx]
      fold_means[is.nan(fold_means) | is.infinite(fold_means)] <- NA_real_

      if (estimation_method == "crossnobis") {
        fold_means[present_idx, ] <- fold_means[present_idx, , drop = FALSE] %*%
          whitening_matrix_W
      }
    }

    non_empty <- rowSums(!is.na(fold_means)) > 0
    if (any(non_empty)) {
      U_hat_sl_cum[non_empty, ] <- U_hat_sl_cum[non_empty, , drop = FALSE] +
        fold_means[non_empty, , drop = FALSE]
      cond_fold_counts[non_empty] <- cond_fold_counts[non_empty] + 1L
    }

    # Store per-fold estimates if requested
    if (return_folds) {
      U_folds_array[, , i] <- fold_means
    }
  }

  # --- Finalise cross-validated means by dividing cumulative sums by counts ---
  divisors <- ifelse(cond_fold_counts > 0, cond_fold_counts, NA_real_)
  U_hat_sl <- sweep(U_hat_sl_cum, 1, divisors, FUN = "/")

  ## ---- Optional post-processing: L2 row normalisation --------------------
  if (estimation_method == "L2_norm") {

    # 1. Compute Euclidean norm of each row (condition pattern)
    row_norms <- sqrt(rowSums(U_hat_sl^2, na.rm = TRUE))

    # 2. Guard against divide-by-zero (all-zero patterns)
    zero_idx <- which(row_norms < 1e-10)
    if(length(zero_idx) > 0){
      warning(sprintf("L2_norm: %d condition patterns had near-zero energy and were left unscaled.",
                      length(zero_idx)), call. = FALSE)
      row_norms[zero_idx] <- 1  # effectively no scaling for those rows
    }

    # 3. Apply scaling (row-wise division)
    U_hat_sl <- sweep(U_hat_sl, 1, row_norms, FUN = "/")
    # If returning folds, L2 norm is applied to the final mean, not per-fold estimates here.
    # The Crossnobis plan implies whitening (if any) happens to per-fold means before distance calc.
    # L2 norm as an estimation_method is about the final U_hat_sl.
  }

  # After U_hat_sl is potentially L2-normalized (if method was "L2_norm"),
  # or is the direct average (if method was "average" or "crossnobis" at this stage before this block).
  if (estimation_method == "average") {
    # For "average" method, if a condition was not present in all folds, 
    # its final mean estimate should be NA.
    if (n_folds > 0) { # Avoid issues if n_folds is 0 for some reason
      for (k_idx in seq_len(n_conditions)) {
        if (cond_fold_counts[k_idx] < n_folds) {
          U_hat_sl[k_idx, ] <- NA_real_
        }
      }
    }
  }

  if (return_folds) {
    return(list(mean_estimate = U_hat_sl, fold_estimates = U_folds_array))
  } else {
    return(U_hat_sl)
  }
}

#' Process a single cross-validation fold
#'
#' Internal helper used by `compute_crossvalidated_means_sl`.
#' It computes per-condition means for one training fold, applies
#' optional whitening, and reindexes the result so that all
#' conditions are represented.
#'
#' @noRd
#' @keywords internal
process_single_fold <- function(train_data_fold,
                               train_labels_fold,
                               unique_conditions,
                               n_conditions,
                               n_voxels,
                               estimation_method,
                               whitening_matrix_W,
                               sl_colnames) {

  fold_means_mat <- NULL
  if (nrow(train_data_fold) > 0 && length(train_labels_fold) > 0) {
    train_labels_factor <- factor(train_labels_fold, levels = unique_conditions)
    sums_by_cond <- rowsum(train_data_fold, group = train_labels_factor,
                           reorder = TRUE, na.rm = TRUE)
    counts_by_cond <- table(train_labels_factor)
    ordered_levels_in_sum <- rownames(sums_by_cond)
    counts_for_division <- counts_by_cond[ordered_levels_in_sum]
    means_by_cond <- sums_by_cond / as.numeric(counts_for_division)
    means_by_cond[is.nan(means_by_cond) | is.infinite(means_by_cond)] <- NA_real_

    fold_means_mat <- means_by_cond
    if (!is.null(colnames(train_data_fold)) &&
        ncol(fold_means_mat) == ncol(train_data_fold)) {
      colnames(fold_means_mat) <- colnames(train_data_fold)
    } else if (!is.null(sl_colnames) &&
               ncol(fold_means_mat) == length(sl_colnames)) {
      colnames(fold_means_mat) <- sl_colnames
    }
  } else {
    fold_means_mat <- matrix(NA_real_, nrow = 0, ncol = n_voxels,
                             dimnames = list(NULL, sl_colnames))
  }

  if (estimation_method == "crossnobis") {
    if (nrow(fold_means_mat) > 0) {
      fold_means_mat_processed <- fold_means_mat %*% whitening_matrix_W
      if (!is.null(colnames(fold_means_mat))) {
        colnames(fold_means_mat_processed) <- colnames(fold_means_mat)
      }
    } else {
      fold_means_mat_processed <- fold_means_mat
    }
  } else {
    fold_means_mat_processed <- fold_means_mat
  }

  full_fold_means <- matrix(NA_real_, nrow = n_conditions, ncol = n_voxels,
                            dimnames = list(unique_conditions, sl_colnames))
  present_conditions <- rownames(fold_means_mat_processed)
  valid_present_conditions <-
    present_conditions[present_conditions %in% rownames(full_fold_means)]
  if (length(valid_present_conditions) > 0) {
    full_fold_means[valid_present_conditions, ] <-
      fold_means_mat_processed[valid_present_conditions, , drop = FALSE]
  }

  full_fold_means
}
