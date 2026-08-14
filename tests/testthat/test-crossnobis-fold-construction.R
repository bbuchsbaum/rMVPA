test_that("crossnobis uses independent blocked partitions", {
  n_folds <- 12L
  n_voxels <- n_folds
  block <- rep(seq_len(n_folds), each = 2L)
  condition <- factor(rep(c("A", "B"), times = n_folds),
                      levels = c("A", "B"))

  # Each partition has an orthogonal condition-difference pattern. Independent
  # cross-partition products are therefore exactly zero, whereas leave-one-out
  # training means share 10 of 12 patterns and yield a positive distance.
  sl_data <- matrix(0, nrow = length(block), ncol = n_voxels,
                    dimnames = list(NULL, paste0("V", seq_len(n_voxels))))
  for (fold in seq_len(n_folds)) {
    sl_data[2L * fold - 1L, fold] <- 0.5
    sl_data[2L * fold, fold] <- -0.5
  }

  mvpa_des <- structure(
    list(
      Y = condition,
      design_matrix = matrix(0, nrow = nrow(sl_data), ncol = 1L)
    ),
    class = c("mvpa_design", "list")
  )
  cv_spec <- blocked_cross_validation(block)

  result <- compute_crossvalidated_means_sl(
    sl_data,
    mvpa_des,
    cv_spec,
    estimation_method = "crossnobis",
    whitening_matrix_W = diag(n_voxels),
    return_folds = TRUE
  )

  expected_folds <- array(
    NA_real_,
    dim = c(2L, n_voxels, n_folds),
    dimnames = list(levels(condition), colnames(sl_data),
                    paste0("Fold", seq_len(n_folds)))
  )
  for (fold in seq_len(n_folds)) {
    expected_folds[, , fold] <- sl_data[block == fold, , drop = FALSE]
  }

  expect_equal(result$fold_estimates, expected_folds, tolerance = 1e-15)

  observed_distance <- rMVPA:::compute_crossnobis_distances_sl(
    result$fold_estimates,
    P_voxels = n_voxels
  )
  expect_equal(unname(observed_distance), 0, tolerance = 1e-15)

  overlapping_folds <- array(
    NA_real_,
    dim = dim(expected_folds),
    dimnames = dimnames(expected_folds)
  )
  for (fold in seq_len(n_folds)) {
    train <- block != fold
    overlapping_folds["A", , fold] <- colMeans(
      sl_data[train & condition == "A", , drop = FALSE]
    )
    overlapping_folds["B", , fold] <- colMeans(
      sl_data[train & condition == "B", , drop = FALSE]
    )
  }

  overlapping_distance <- rMVPA:::compute_crossnobis_distances_sl(
    overlapping_folds,
    P_voxels = n_voxels
  )
  expected_bias <- (n_folds - 2) / ((n_folds - 1)^2 * n_voxels)
  expect_equal(unname(overlapping_distance), expected_bias, tolerance = 1e-15)
  expect_gt(unname(overlapping_distance), 0)
})
