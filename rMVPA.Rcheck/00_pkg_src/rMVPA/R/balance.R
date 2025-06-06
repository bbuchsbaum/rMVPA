#' Balance Cross-Validation Partitions
#'
#' Modifies a cross-validation partitioning scheme to ensure that each
#' target class has an equal number of samples within each training fold,
#' and optionally within each test fold, using either sub-sampling or
#' oversampling.
#'
#' @param obj A `cross_validation` object (e.g., from
#'   `blocked_cross_validation`, `kfold_cross_validation`).
#' @param design An `mvpa_design` object containing the target labels
#'   (`.sa.targets`) corresponding to the original dataset samples.
#' @param method Character string specifying the balancing method:
#'   - `"subsample"` (default): Down-samples majority classes to match the
#'     size of the smallest class (sampling without replacement).
#'   - `"oversample"`: Up-samples minority classes to match the size of the
#'     largest class (sampling with replacement).
#' @param ... Additional arguments passed to internal balancing functions:
#'   \describe{
#'     \item{balance_test_set}{Logical. If `TRUE` (default), balance the test
#'       sets in each fold as well using the specified `method`. **Note:**
#'       Oversampling the test set is generally not recommended as it can
#'       lead to misleading performance estimates. A warning will be issued if
#'       `balance_test_set=TRUE` and `method="oversample"`.}
#'     \item{seed}{An optional integer seed for the random number generator for
#'       reproducible sampling. If `NULL` (default), the result varies.}
#'   }
#'
#' @return A `custom_cross_validation` object where the sample indices in
#'   `.train_indices` (and optionally `.test_indices`) for each fold have
#'   been resampled to ensure balanced representation of target classes.
#'
#' @details
#'
#' **Sub-sampling (`method="subsample"`)**:
#' 1. Identifies the target class with the minimum number of samples
#'    (`min_count`) within the set (train or test).
#' 2. For *every* target class within that set, it randomly selects exactly
#'    `min_count` samples *without replacement*.
#' 3. Discards samples from majority classes.
#'
#' **Oversampling (`method="oversample"`)**:
#' 1. Identifies the target class with the maximum number of samples
#'    (`max_count`) within the set (train or test).
#' 2. For *every* target class within that set, it randomly selects exactly
#'    `max_count` samples *with replacement*.
#' 3. Duplicates samples from minority classes.
#'
#' Balancing guarantees that after the process, each target class appears
#' equally often within each balanced training (and optionally testing) set.
#' This is useful for preventing classifiers from being biased towards
#' majority classes.
#'
#' The output is always a `custom_cross_validation` object because the
#' balancing process defines specific, explicit sets of indices for each
#' fold, which may no longer strictly adhere to the original blocking or
#' k-fold structure.
#'
#' @examples
#' # Create an imbalanced dataset design (more class 'b')
#' design_df <- data.frame(condition = factor(rep(c("a", "b", "b"), 20)),
#'                        block = rep(1:6, each = 10))
#' des <- mvpa_design(design_df, y_train = ~ condition, block_var = ~ block)
#'
#' # Create standard blocked partitions (likely unbalanced)
#' cval_unbalanced <- blocked_cross_validation(des$block_var)
#' print("Unbalanced Counts (Example Fold 1 Train):")
#' print(table(des$y_train[unlist(crossval_samples(cval_unbalanced,
#'           design_df, des$y_train)$train[[1]]$idx)]))
#'
#' # Balance using sub-sampling (default)
#' cval_sub <- balance_partitions(cval_unbalanced, des, seed = 1)
#' print(cval_sub)
#' print("Subsample Balanced Counts (Example Fold 1 Train):")
#' print(table(crossval_samples(cval_sub, design_df, des$y_train)$ytrain[[1]]))
#'
#' # Balance using over-sampling
#' cval_over <- balance_partitions(cval_unbalanced, des, method = "oversample", seed = 2)
#' print(cval_over)
#' print("Oversample Balanced Counts (Example Fold 1 Train):")
#' print(table(crossval_samples(cval_over, design_df, des$y_train)$ytrain[[1]]))
#'
#' @export
#' @family cross_validation
balance_partitions <- function(obj, design, ...) {
  UseMethod("balance_partitions")
}


#' @rdname balance_partitions
#' @export
balance_partitions.default <- function(obj, design, method = "subsample", ...) {
  warning("Balancing not explicitly implemented for class '",
          class(obj)[1], "'. Attempting default balancing using method='", method, "'.", call. = FALSE)
  perform_balancing_on_cv_obj(obj, design, method = method, ...)
}

#' @rdname balance_partitions
#' @export
balance_partitions.blocked_cross_validation <- function(obj, design, method = "subsample", ...) {
  perform_balancing_on_cv_obj(obj, design, method = method, ...)
}

#' @rdname balance_partitions
#' @export
balance_partitions.kfold_cross_validation <- function(obj, design, method = "subsample", ...) {
  perform_balancing_on_cv_obj(obj, design, method = method, ...)
}

#' @rdname balance_partitions
#' @export
balance_partitions.twofold_blocked_cross_validation <- function(obj, design, method = "subsample", ...) {
  perform_balancing_on_cv_obj(obj, design, method = method, ...)
}

#' @rdname balance_partitions
#' @export
balance_partitions.bootstrap_blocked_cross_validation <- function(obj, design, method = "subsample", ...) {
  warning("Balancing bootstrap partitions might have unintended consequences. ",
          "Bootstrap already involves resampling. Proceeding with method='", method, "'.", call. = FALSE)
  perform_balancing_on_cv_obj(obj, design, method = method, ...)
}

#' @rdname balance_partitions
#' @export
balance_partitions.sequential_blocked_cross_validation <- function(obj, design, method = "subsample", ...) {
  perform_balancing_on_cv_obj(obj, design, method = method, ...)
}

#' @rdname balance_partitions
#' @export
balance_partitions.custom_cross_validation <- function(obj, design, method = "subsample", ...) {
  # Custom CV might already be balanced, but we allow re-balancing
  perform_balancing_on_cv_obj(obj, design, method = method, ...)
}

# Internal helper function to perform the balancing logic
# Takes a cross_validation object and design, returns a custom_cross_validation obj
#' @keywords internal
#' @noRd
perform_balancing_on_cv_obj <- function(obj, design, method = "subsample",
                                        balance_test_set = TRUE, seed = NULL, ...) {

  # Validate method argument
  valid_methods <- c("subsample", "oversample")
  if (!method %in% valid_methods) {
    stop("`method` must be one of: ", paste(valid_methods, collapse=", "), call. = FALSE)
  }

  if (!inherits(design, "mvpa_design")) {
    stop("`design` argument must be an 'mvpa_design' object.")
  }

  # Issue warning if oversampling the test set
  if (balance_test_set && method == "oversample") {
      warning("Oversampling the test set ('balance_test_set = TRUE', ",
              "method = 'oversample') is generally not recommended and ",
              "may lead to inflated performance metrics.", call. = FALSE)
  }


  targets <- design$y_train # Use training targets for reference

  # Get original sample indices for each fold using the appropriate method
  n_obs_total <- length(targets)
  dummy_data <- data.frame(dummy_col = seq_len(n_obs_total))
  initial_folds_tibble <- crossval_samples(obj, dummy_data, targets)

  balanced_sample_set <- vector("list", nrow(initial_folds_tibble))

  # Set seed if provided
  if (!is.null(seed)) {
    if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) runif(1)
    old_seed <- .GlobalEnv$.Random.seed
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }

  # Iterate through each fold provided by crossval_samples
  for (i in 1:nrow(initial_folds_tibble)) {
    train_indices_orig <- initial_folds_tibble$train[[i]]$idx
    test_indices_orig  <- initial_folds_tibble$test[[i]]$idx

    # Balance training set
    balanced_train_indices <- balance_single_fold(train_indices_orig, targets, method)

    # Balance test set if requested
    balanced_test_indices <- if (balance_test_set) {
      balance_single_fold(test_indices_orig, targets, method)
    } else {
      test_indices_orig # Keep original test indices
    }

    # Store the balanced train/test indices for this fold
    balanced_sample_set[[i]] <- list(train = balanced_train_indices,
                                     test = balanced_test_indices)
  }

  # Return a new custom_cross_validation object with the balanced indices
  custom_cross_validation(balanced_sample_set)
}


# Internal helper function to balance indices for a single fold
# Now includes 'method' argument
#' @keywords internal
#' @noRd
balance_single_fold <- function(fold_indices, all_targets, method) {
  if (length(fold_indices) == 0) {
    return(integer(0)) # Return empty if fold is empty
  }

  fold_targets <- all_targets[fold_indices]
  target_counts <- table(fold_targets)

  if (length(target_counts) <= 1) {
      # Only one class present or empty fold, cannot balance further
      warning("Fold contains only one class or is empty; cannot balance.", call. = FALSE)
      return(fold_indices) # Return original indices for this fold
  }

  # Check for classes with zero samples in this fold (should ideally not happen
  # if partitions are valid, but handle defensively)
  if (any(target_counts == 0)) {
      warning("Fold contains a class with zero samples; resulting balanced fold might be incomplete or invalid.", call.=FALSE)
      # Proceed, but the target size calculation might be affected
      target_counts <- target_counts[target_counts > 0] # Exclude zero-count classes
      if (length(target_counts) <= 1) {
          return(fold_indices) # Cannot balance if only one class remains
      }
  }


  balanced_indices <- integer(0)

  # Get unique targets present in this specific fold
  unique_fold_targets <- names(target_counts)

  if (method == "subsample") {
    # --- Sub-sampling logic ---
    target_size <- min(target_counts)

    for (target_level in unique_fold_targets) {
      indices_this_target_in_fold <- fold_indices[which(fold_targets == target_level)]
      n_available <- length(indices_this_target_in_fold)
      # sample size is target_size (min_count)
      sample_size <- min(target_size, n_available) # Defensive check

      if (n_available > 0) {
          sampled_indices <- if (n_available == 1 && sample_size == 1) {
                                 indices_this_target_in_fold
                             } else {
                                 sample(indices_this_target_in_fold, size = sample_size, replace = FALSE)
                             }
          balanced_indices <- c(balanced_indices, sampled_indices)
      }
    }

  } else if (method == "oversample") {
    # --- Oversampling logic ---
    target_size <- max(target_counts)

    for (target_level in unique_fold_targets) {
      indices_this_target_in_fold <- fold_indices[which(fold_targets == target_level)]
      n_available <- length(indices_this_target_in_fold)

      if (n_available > 0) {
          # sample size is target_size (max_count)
          # Sample *with replacement*
          sampled_indices <- sample(indices_this_target_in_fold, size = target_size, replace = TRUE)
          balanced_indices <- c(balanced_indices, sampled_indices)
      } # If n_available is 0, this class cannot be oversampled
    }

  } else {
      # Should have been caught earlier, but double-check
      stop("Internal error: invalid method '", method, "' in balance_single_fold.")
  }

  # Return the sorted balanced indices for consistency
  sort(balanced_indices)
}