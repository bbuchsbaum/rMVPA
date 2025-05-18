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
#'   This matrix should be the whitening transformation (e.g., Σ_noise^(-1/2)) derived from GLM residuals.
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
  if (!is.matrix(sl_data) || !is.numeric(sl_data)) {
    rlang::abort("`sl_data` must be a numeric matrix (samples x features).")
  }
  if (!inherits(mvpa_design, "mvpa_design")) {
    rlang::abort("`mvpa_design` must be an object of class 'mvpa_design'.")
  }
  
  if (!inherits(cv_spec, "cross_validation")) {
    rlang::abort("`cv_spec` must be an object inheriting from class 'cross_validation'.")
  }
  if (nrow(sl_data) != nrow(mvpa_design$design_matrix)) {
     # Assuming design_matrix rows correspond to samples
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
  # Iterate folds and accumulate means online (memory safe)
  for (i in seq_len(n_folds)) {

    train_indices <- train_indices(cv_spec, i)

    if (length(train_indices) == 0) {
      warning(paste("Fold", i, "has no training samples. Skipping."))
      next
    }

    train_data_fold <- sl_data[train_indices, , drop = FALSE]
    train_labels_fold <- condition_labels[train_indices]

    fold_means_mat <- NULL
    # Ensure there are rows and labels to process
    if (nrow(train_data_fold) > 0 && length(train_labels_fold) > 0) {
        # Convert train_labels_fold to a factor to ensure consistent ordering for rowsum and table
        # Drop unused levels if any, to prevent issues with table counts vs rowsum rows
        train_labels_factor <- factor(train_labels_fold, levels = unique_conditions)
        
        # Calculate sums per condition using rowsum
        # rowsum requires the group to be a factor for predictable row order in output
        # It sorts by levels of the factor.
        sums_by_cond <- rowsum(train_data_fold, group = train_labels_factor, reorder = TRUE, na.rm = TRUE)
        
        # Calculate counts per condition
        counts_by_cond <- table(train_labels_factor)
        
        # Align counts with the order of sums_by_cond (which is sorted by factor levels)
        # This ensures correct element-wise division.
        # If a level had sum 0 (all NAs) but count > 0, or sum > 0 but count 0 (should not happen here if table from same factor)
        # we need to be careful.
        # Rownames of sums_by_cond are the factor levels that had data.
        ordered_levels_in_sum <- rownames(sums_by_cond)
        counts_for_division <- counts_by_cond[ordered_levels_in_sum]
        
        # Calculate means: sums / counts
        # Avoid division by zero if a count is zero (should not happen if level in sums_by_cond)
        means_by_cond <- sums_by_cond / as.numeric(counts_for_division)
        # Handle cases where a sum might have been NA (if all inputs were NA) or count was 0 (resulting in NaN/Inf)
        means_by_cond[is.nan(means_by_cond) | is.infinite(means_by_cond)] <- NA_real_
        
        fold_means_mat <- means_by_cond
        # rownames are already correctly set by rowsum to be the condition labels (factor levels)
        # Ensure colnames are preserved from original data if possible (rowsum might strip them)
        if (!is.null(colnames(train_data_fold)) && ncol(fold_means_mat) == ncol(train_data_fold)) {
            colnames(fold_means_mat) <- colnames(train_data_fold)
        } else if (!is.null(colnames(sl_data)) && ncol(fold_means_mat) == ncol(sl_data)) {
            colnames(fold_means_mat) <- colnames(sl_data) # Fallback to sl_data colnames
        }

    } else {
        # If no training data/labels for this fold, create an empty matrix with correct colnames
        # This ensures consistency for later re-indexing into full_fold_means
        num_unique_conds_total <- length(unique_conditions) # from outer scope
        k_fold_levels <- character(0)
        if (length(train_labels_fold) > 0) k_fold_levels <- levels(factor(train_labels_fold)) 
        # If no data, it means no conditions from this fold, so 0 rows.
        # The re-indexing to full_fold_means will handle this by having all NAs for this fold.
        fold_means_mat <- matrix(NA_real_, 
                                 nrow = 0, 
                                 ncol = ncol(sl_data), 
                                 dimnames = list(NULL, colnames(sl_data)))
    }

    # Apply whitening if method is crossnobis
    if (estimation_method == "crossnobis") {
        # We already checked that whitening_matrix_W is valid and has correct dimensions (V_sl x V_sl)
        # fold_means_mat is K_fold x V_sl (where K_fold is conditions present in this fold)
        if (nrow(fold_means_mat) > 0) { # Ensure there are rows to whiten
             fold_means_mat_processed <- fold_means_mat %*% whitening_matrix_W
             # Ensure column names are preserved if they exist on W (though not strictly necessary here)
             if (!is.null(colnames(whitening_matrix_W)) && ncol(fold_means_mat_processed) == ncol(whitening_matrix_W)){
                # this might not be right if W's colnames are not feature names
             } else if (!is.null(colnames(fold_means_mat))){
                colnames(fold_means_mat_processed) <- colnames(fold_means_mat)
             }
        } else {
             fold_means_mat_processed <- fold_means_mat # No rows to process
        }
    } else {
        fold_means_mat_processed <- fold_means_mat # Other methods don't whiten here
    }

    # Handle cases where a condition might be missing in a fold's training set
    # Reindex to ensure all conditions are present, filling missing with NA
    full_fold_means <- matrix(NA_real_,
                              nrow = n_conditions,
                              ncol = n_voxels,
                              dimnames = list(unique_conditions, colnames(sl_data)))
                              
    present_conditions <- rownames(fold_means_mat_processed) # Use processed matrix
    # Ensure present_conditions are valid before trying to subset
    valid_present_conditions <- present_conditions[present_conditions %in% rownames(full_fold_means)]
    if (length(valid_present_conditions) > 0) {
        full_fold_means[valid_present_conditions, ] <- fold_means_mat_processed[valid_present_conditions, , drop = FALSE]
    }
    
    # Update cumulative matrix and counts
    # Only add to cumulative sum if the full_fold_means for that condition in that fold is not all NA
    # This handles cases where a condition might be entirely missing from a fold, 
    # and its full_fold_means row would be all NAs.
    for (k_idx in seq_len(n_conditions)) {
        if (any(!is.na(full_fold_means[k_idx, ]))) { # If at least one non-NA value in the row
            # Add to cumulative, replacing NAs in cumulative with the new values if cum was 0
            # or adding to existing values. This needs care if U_hat_sl_cum starts at 0.
            # If U_hat_sl_cum[k_idx,] is all 0s from init, and full_fold_means[k_idx,] has NAs,
            # 0 + NA = NA. This is correct.
            U_hat_sl_cum[k_idx, ] <- U_hat_sl_cum[k_idx, ] + full_fold_means[k_idx, ]
            cond_fold_counts[k_idx] <- cond_fold_counts[k_idx] + 1
        }
    }

    # Store per-fold estimates if requested
    if (return_folds) {
      U_folds_array[, , i] <- full_fold_means
    }
  }

  # --- Finalise cross-validated means by dividing cumulative sums by counts ---
  safe_div <- function(sum_vec, count) {
       ifelse(count > 0, sum_vec / count, NA_real_)
  }

  U_hat_sl <- U_hat_sl_cum
  for (k in seq_len(n_conditions)) {
       U_hat_sl[k, ] <- safe_div(U_hat_sl_cum[k, ], cond_fold_counts[k])
  }

  # Rows where cond_fold_counts == 0 become all NA (indicating condition never in training)
  zero_rows <- which(cond_fold_counts == 0)
  if (length(zero_rows) > 0) {
      U_hat_sl[zero_rows, ] <- NA_real_
  }

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