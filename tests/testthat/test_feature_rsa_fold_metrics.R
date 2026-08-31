library(rMVPA)

fold_metric_reference <- function(predicted, observed, fold_id) {
  cross_cor <- suppressWarnings(stats::cor(
    t(predicted),
    t(observed),
    use = "pairwise.complete.obs"
  ))
  diagonal <- diag(cross_cor)
  ranks <- rep(NA_real_, nrow(cross_cor))
  off_diagonal <- numeric(0)

  for (i in seq_len(nrow(cross_cor))) {
    candidates <- which(fold_id == fold_id[[i]])
    row_values <- cross_cor[i, candidates]
    correct <- match(i, candidates)
    denominator <- sum(is.finite(row_values)) - 1L
    if (denominator > 0L && is.finite(row_values[[correct]])) {
      ranks[[i]] <-
        (sum(row_values <= row_values[[correct]], na.rm = TRUE) - 1L) /
        denominator
    }
    keep <- is.finite(row_values)
    keep[[correct]] <- FALSE
    off_diagonal <- c(off_diagonal, row_values[keep])
  }

  pattern_correlation <- mean(diagonal, na.rm = TRUE)
  c(
    correlation = pattern_correlation,
    discrimination = pattern_correlation - mean(off_diagonal),
    rank = mean(ranks, na.rm = TRUE)
  )
}

fold_rdm_reference <- function(x, fold_id) {
  correlation <- suppressWarnings(stats::cor(
    t(x),
    use = "pairwise.complete.obs"
  ))
  values <- as.numeric((1 - correlation)[lower.tri(correlation)])
  within_fold <- outer(fold_id, fold_id, FUN = "==")[lower.tri(correlation)]
  values[!within_fold] <- NA_real_
  values
}

test_that("fold-aware pattern metrics match a dense grouped oracle", {
  set.seed(7701)
  observed <- matrix(rnorm(18 * 24), 18, 24)
  predicted <- observed * 0.45 + matrix(rnorm(18 * 24), 18, 24)
  fold_id <- rep(c("film-c", "film-a", "film-b"), each = 6L)

  expected <- fold_metric_reference(predicted, observed, fold_id)
  got <- rMVPA:::.feature_rsa_pattern_metrics_blockwise(
    predicted,
    observed,
    fold_id = fold_id
  )

  expect_equal(got, expected, tolerance = 2e-12)

  withr::local_options(rMVPA.feature_rsa_metric_block_rows = 1L)
  one_row_blocks <- rMVPA:::.feature_rsa_pattern_metrics_blockwise(
    predicted,
    observed,
    fold_id = fold_id
  )
  expect_equal(one_row_blocks, expected, tolerance = 2e-12)

  order <- sample(seq_len(nrow(observed)))
  reordered <- rMVPA:::.feature_rsa_pattern_metrics_blockwise(
    predicted[order, , drop = FALSE],
    observed[order, , drop = FALSE],
    fold_id = paste0("renamed-", match(fold_id[order], unique(fold_id[order])))
  )
  expect_equal(reordered, expected, tolerance = 2e-12)
})

test_that("intercept-only OOF predictions have chance rank within held-out folds", {
  set.seed(7702)
  n <- 60L
  n_voxels <- 100L
  fold_id <- rep(seq_len(3L), each = n / 3L)
  observed <- matrix(rnorm(n * n_voxels), n, n_voxels)
  predicted <- matrix(NA_real_, n, n_voxels)

  for (fold in unique(fold_id)) {
    test <- which(fold_id == fold)
    predicted[test, ] <- matrix(
      colMeans(observed[-test, , drop = FALSE]),
      nrow = length(test),
      ncol = n_voxels,
      byrow = TRUE
    )
  }

  pooled <- evaluate_model.feature_rsa_model(
    object = NULL,
    predicted = predicted,
    observed = observed
  )
  heldout <- evaluate_model.feature_rsa_model(
    object = NULL,
    predicted = predicted,
    observed = observed,
    fold_id = fold_id
  )

  expect_lt(pooled$pattern_rank_percentile, 0.45)
  expect_equal(heldout$pattern_rank_percentile, 0.5, tolerance = 2e-12)
  expect_equal(heldout$pattern_discrimination, 0, tolerance = 2e-12)
  expect_equal(
    heldout$pattern_correlation,
    pooled$pattern_correlation,
    tolerance = 2e-12
  )
})

test_that("fold-aware RDM vectors retain shape and mask cross-fold pairs", {
  set.seed(7703)
  observed <- matrix(rnorm(15 * 20), 15, 20)
  predicted <- observed + matrix(rnorm(15 * 20, sd = 0.35), 15, 20)
  fold_id <- rep(letters[1:3], each = 5L)
  expected_predicted <- fold_rdm_reference(predicted, fold_id)
  expected_observed <- fold_rdm_reference(observed, fold_id)

  result <- evaluate_model.feature_rsa_model(
    object = NULL,
    predicted = predicted,
    observed = observed,
    fold_id = fold_id,
    compute_rdm_vectors = TRUE
  )

  expect_equal(
    result$predicted_rdm_vec,
    expected_predicted,
    tolerance = 2e-12
  )
  expect_equal(
    result$observed_rdm_vec,
    expected_observed,
    tolerance = 2e-12
  )
  expect_length(result$predicted_rdm_vec, choose(nrow(observed), 2L))
  expect_equal(sum(is.finite(result$predicted_rdm_vec)), 3L * choose(5L, 2L))
  expect_equal(
    result$rdm_correlation,
    suppressWarnings(stats::cor(
      expected_predicted,
      expected_observed,
      method = "spearman",
      use = "complete.obs"
    )),
    tolerance = 2e-12
  )
})

test_that("fold-aware geometry validates groups and handles singleton folds", {
  set.seed(7704)
  observed <- matrix(rnorm(8 * 12), 8, 12)
  predicted <- observed + matrix(rnorm(8 * 12), 8, 12)

  expect_error(
    evaluate_model.feature_rsa_model(
      NULL,
      predicted,
      observed,
      fold_id = 1:7
    ),
    "fold_id.*8"
  )
  expect_error(
    evaluate_model.feature_rsa_model(
      NULL,
      predicted,
      observed,
      fold_id = c(1:7, NA)
    ),
    "fold_id.*missing"
  )

  singleton <- evaluate_model.feature_rsa_model(
    NULL,
    predicted,
    observed,
    fold_id = seq_len(nrow(observed))
  )
  expect_true(is.finite(singleton$pattern_correlation))
  expect_true(is.na(singleton$pattern_discrimination))
  expect_true(is.na(singleton$pattern_rank_percentile))
  expect_true(is.na(singleton$rdm_correlation))
})

test_that("feature RSA variance kernels are stable for constant and shifted data", {
  set.seed(7705)
  base <- matrix(rnorm(40 * 13), 40, 13)
  repeated <- matrix(colMeans(base), 25, ncol(base), byrow = TRUE)

  expect_identical(
    rMVPA:::.feature_rsa_col_sds(repeated),
    rep(0, ncol(repeated))
  )
  expect_identical(
    rMVPA:::.feature_rsa_row_sds(matrix(7.25, 11, 9)),
    rep(0, 11L)
  )

  shifted_columns <- sweep(base, 2L, seq(1e10, 1.3e10, length.out = ncol(base)), "+")
  shifted_rows <- sweep(base, 1L, seq(1e10, 1.3e10, length.out = nrow(base)), "+")
  expect_equal(
    rMVPA:::.feature_rsa_col_sds(shifted_columns),
    apply(shifted_columns, 2L, stats::sd),
    tolerance = 1e-12,
    scale = 1
  )
  expect_equal(
    rMVPA:::.feature_rsa_row_sds(shifted_rows),
    apply(shifted_rows, 1L, stats::sd),
    tolerance = 1e-12,
    scale = 1
  )
})

test_that("fold permutation indices never cross outer-fold boundaries", {
  fold_id <- c(
    rep(c("run-3", "run-1", "run-2"), c(3L, 5L, 4L)),
    "singleton"
  )
  set.seed(7706)
  for (iteration in seq_len(50L)) {
    index <- rMVPA:::.feature_rsa_permutation_index(fold_id)
    expect_setequal(index, seq_along(fold_id))
    expect_identical(fold_id[index], fold_id)
    expect_identical(index[[length(index)]], length(index))
  }
})
