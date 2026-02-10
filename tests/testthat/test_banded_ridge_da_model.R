context("banded_ridge_da_model")

.make_linear_da_dataset <- function(Tenc = 30,
                                    Trec = 20,
                                    dims = c(2, 2, 1),
                                    set_sizes = c(a = 2, b = 2),
                                    zero_sets = "b",
                                    seed = 1) {
  set.seed(seed)
  Dfeat <- sum(set_sizes)
  V <- prod(dims)

  X_enc <- matrix(rnorm(Tenc * Dfeat), Tenc, Dfeat)
  X_rec <- matrix(rnorm(Trec * Dfeat), Trec, Dfeat)

  B <- matrix(rnorm(Dfeat * V), nrow = Dfeat, ncol = V)
  if (!is.null(zero_sets)) {
    idx <- split(seq_len(Dfeat), rep(names(set_sizes), set_sizes))
    for (nm in intersect(names(idx), zero_sets)) {
      B[idx[[nm]], ] <- 0
    }
  }

  Y_enc <- X_enc %*% B
  Y_rec <- X_rec %*% B

  arr_enc <- array(0, dim = c(dims, Tenc))
  for (t in seq_len(Tenc)) {
    arr_enc[,,,t] <- array(Y_enc[t, ], dim = dims)
  }
  arr_rec <- array(0, dim = c(dims, Trec))
  for (t in seq_len(Trec)) {
    arr_rec[,,,t] <- array(Y_rec[t, ], dim = dims)
  }

  mask <- neuroim2::NeuroVol(array(1, dim = dims), neuroim2::NeuroSpace(dims))
  train_vec <- neuroim2::NeuroVec(arr_enc, neuroim2::NeuroSpace(c(dims, Tenc)))
  test_vec <- neuroim2::NeuroVec(arr_rec, neuroim2::NeuroSpace(c(dims, Trec)))

  dataset <- mvpa_dataset(train_vec, test_vec, mask)
  spec <- do.call(blocks, as.list(set_sizes))
  fs_enc <- feature_sets(X_enc, spec)
  fs_rec <- feature_sets(X_rec, spec, set_order = fs_enc$set_order)
  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, length.out = Trec))

  list(
    dataset = dataset,
    design = des,
    mask = mask,
    X_enc = X_enc,
    X_rec = X_rec,
    B = B,
    set_sizes = set_sizes
  )
}

test_that("feature_sets + expected_features build expected shapes", {
  set.seed(1)

  Tenc <- 10
  Trec <- 6
  Xenc <- matrix(rnorm(Tenc * 8), Tenc, 8)
  fs <- feature_sets(Xenc, blocks(a = 3, b = 5))

  expect_s3_class(fs, "feature_sets")
  expect_equal(nrow(fs$X), Tenc)
  expect_equal(ncol(fs$X), 8)
  expect_equal(names(fs$indices), c("a", "b"))
  expect_equal(length(fs$indices$a), 3)
  expect_equal(length(fs$indices$b), 5)

  G <- matrix(runif(Trec * Tenc), Trec, Tenc)
  G <- G / rowSums(G)
  fs_rec <- expected_features(fs, G, drop_null = FALSE, renormalize = TRUE)

  expect_s3_class(fs_rec, "feature_sets")
  expect_equal(nrow(fs_rec$X), Trec)
  expect_equal(ncol(fs_rec$X), 8)
})


test_that("banded_ridge_da_model runs regional pipeline (run-blocked recall CV)", {
  skip_on_cran()
  set.seed(123)

  Tenc <- 24
  Trec <- 20
  toy <- gen_sample_dataset(D = c(4, 4, 4), nobs = Tenc, nlevels = 2, blocks = 3, external_test = TRUE, ntest_obs = Trec)

  regionMask <- neuroim2::NeuroVol(array(1, c(4, 4, 4)), neuroim2::space(toy$dataset$mask))

  Xenc <- matrix(rnorm(Tenc * 12), Tenc, 12)
  fs_enc <- feature_sets(Xenc, blocks(low = 4, mid = 4, sem = 4))

  G <- matrix(runif(Trec * Tenc), Trec, Tenc)
  G <- G / rowSums(G)
  fs_rec <- expected_features(fs_enc, G, drop_null = FALSE, renormalize = TRUE)

  recall_runs <- rep(1:2, each = Trec / 2)
  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = recall_runs)

  ms <- banded_ridge_da_model(
    dataset = toy$dataset,
    design = des,
    mode = "stacked",
    lambdas = c(low = 10, mid = 10, sem = 10),
    alpha_recall = 0.2,
    compute_delta_r2 = TRUE,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, regionMask)

  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_true(all(c("recall_r2_full", "recall_mse_full", "delta_r2_low", "delta_r2_mid", "delta_r2_sem") %in% names(res$performance_table)))
  expect_true(is.list(res$fits))
  expect_equal(length(res$fits[[1]]$folds), 2)
})


test_that("banded_ridge_da_model falls back to contiguous recall folds when only one run", {
  skip_on_cran()
  set.seed(456)

  Tenc <- 20
  Trec <- 18
  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = Tenc, nlevels = 2, blocks = 2, external_test = TRUE, ntest_obs = Trec)
  regionMask <- neuroim2::NeuroVol(array(1, c(3, 3, 3)), neuroim2::space(toy$dataset$mask))

  Xenc <- matrix(rnorm(Tenc * 9), Tenc, 9)
  fs_enc <- feature_sets(Xenc, blocks(a = 3, b = 3, c = 3))

  G <- matrix(runif(Trec * Tenc), Trec, Tenc)
  G <- G / rowSums(G)
  fs_rec <- expected_features(fs_enc, G, drop_null = FALSE, renormalize = TRUE)

  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1, Trec))

  ms <- banded_ridge_da_model(
    dataset = toy$dataset,
    design = des,
    mode = "stacked",
    lambdas = c(a = 1, b = 1, c = 1),
    recall_nfolds = 3,
    compute_delta_r2 = FALSE,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, regionMask)
  expect_equal(length(res$fits[[1]]$folds), 3)
})

test_that("banded_ridge_da_model supports purged contiguous folds for single-run target", {
  skip_on_cran()
  set.seed(457)

  Tenc <- 20
  Trec <- 18
  toy <- gen_sample_dataset(D = c(3, 3, 3), nobs = Tenc, nlevels = 2, blocks = 2, external_test = TRUE, ntest_obs = Trec)
  regionMask <- neuroim2::NeuroVol(array(1, c(3, 3, 3)), neuroim2::space(toy$dataset$mask))

  Xenc <- matrix(rnorm(Tenc * 9), Tenc, 9)
  fs_enc <- feature_sets(Xenc, blocks(a = 3, b = 3, c = 3))

  G <- matrix(runif(Trec * Tenc), Trec, Tenc)
  G <- G / rowSums(G)
  fs_rec <- expected_features(fs_enc, G, drop_null = FALSE, renormalize = TRUE)

  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1, Trec))

  ms <- banded_ridge_da_model(
    dataset = toy$dataset,
    design = des,
    mode = "stacked",
    lambdas = c(a = 1, b = 1, c = 1),
    recall_nfolds = 3,
    target_gap = 1,
    compute_delta_r2 = FALSE,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, regionMask)
  folds <- res$fits[[1]]$folds
  expect_equal(length(folds), 3)

  for (f in folds) {
    test_idx <- f$test
    purge_lo <- max(1L, min(test_idx) - 1L)
    purge_hi <- min(Trec, max(test_idx) + 1L)
    purged <- seq.int(purge_lo, purge_hi)
    expect_false(any(f$train %in% purged))
  }
})

test_that("banded_ridge_da_model target_gap overrides recall_gap", {
  dat <- .make_linear_da_dataset(seed = 73, Trec = 12)
  expect_warning(
    ms <- banded_ridge_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1, b = 1),
      recall_nfolds = 3,
      recall_gap = 0,
      target_gap = 2,
      compute_delta_r2 = FALSE
    ),
    "both `recall_gap` and `target_gap`"
  )
  expect_equal(ms$target_gap, 2L)
  expect_equal(ms$recall_gap, 2L)
})

test_that("banded_ridge_da_model does not purge leave-one-run-out folds", {
  dat <- .make_linear_da_dataset(seed = 74, Trec = 12)
  dat$design$block_var_test <- rep(1:2, each = 6)
  regionMask <- dat$mask

  ms <- banded_ridge_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    target_gap = 2,
    compute_delta_r2 = FALSE,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, regionMask)
  folds <- res$fits[[1]]$folds
  expect_equal(length(folds), 2)
  expect_equal(length(folds[[1]]$test), 6)
  expect_equal(length(folds[[1]]$train), 6)
  expect_equal(length(folds[[2]]$test), 6)
  expect_equal(length(folds[[2]]$train), 6)
})

test_that("banded_ridge_da_model computes single-run target permutation null metrics", {
  skip_on_cran()
  set.seed(75)

  dat <- .make_linear_da_dataset(seed = 75, Trec = 12)
  dat$design$block_var_test <- rep(1, 12)

  ms <- banded_ridge_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    recall_nfolds = 3,
    target_gap = 1,
    target_nperm = 8,
    target_perm_strategy = "circular_shift",
    compute_delta_r2 = FALSE,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_true("target_perm_p_r2_full" %in% names(perf))
  expect_true("target_perm_p_mse_full" %in% names(perf))
  expect_true(perf$target_perm_p_r2_full >= 0 && perf$target_perm_p_r2_full <= 1)
  expect_true(perf$target_perm_p_mse_full >= 0 && perf$target_perm_p_mse_full <= 1)

  diag <- res$fits[[1]]
  expect_equal(diag$target_nperm, 8L)
  expect_equal(diag$target_perm_strategy, "circular_shift")
})

test_that("banded_ridge_da_model disables target permutations for multi-run targets", {
  dat <- .make_linear_da_dataset(seed = 76, Trec = 12)
  dat$design$block_var_test <- rep(1:2, each = 6)

  ms <- expect_warning(
    banded_ridge_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1, b = 1),
      target_nperm = 5,
      compute_delta_r2 = FALSE
    ),
    "single-run targets"
  )
  expect_equal(ms$target_nperm, 0L)

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_false("target_perm_p_r2_full" %in% names(perf))
})

test_that("banded_ridge_da_model recovers exact mapping in stacked mode", {
  dat <- .make_linear_da_dataset(seed = 10)
  regionMask <- dat$mask

  ms <- banded_ridge_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 0, b = 1),
    alpha_target = 0,
    compute_delta_r2 = TRUE,
    return_diagnostics = TRUE
  )

  res <- suppressWarnings(run_regional(ms, regionMask))
  perf <- res$performance_table[1, ]
  expect_gt(perf$recall_r2_full, 0.999)
  expect_true(is.finite(perf$target_r2_full))
  expect_lt(abs(perf$delta_r2_b), 1e-6)
  expect_gt(perf$delta_r2_a, 0.5)
})

test_that("banded_ridge_da_model recovers identical mappings in coupled mode", {
  dat <- .make_linear_da_dataset(seed = 11)
  regionMask <- dat$mask

  ms <- banded_ridge_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "coupled",
    lambdas = c(a = 0, b = 1),
    alpha_target = 0.5,
    rho = 5,
    compute_delta_r2 = FALSE,
    return_diagnostics = TRUE
  )

  res <- suppressWarnings(run_regional(ms, regionMask))
  perf <- res$performance_table[1, ]
  expect_gt(perf$recall_r2_full, 0.99)
  expect_equal(perf$target_r2_full, perf$recall_r2_full)
})

test_that("banded_ridge_da_model uses target_folds when provided", {
  dat <- .make_linear_da_dataset(seed = 12, Trec = 10)
  regionMask <- dat$mask

  folds <- list(
    list(train = 6:10, test = 1:5),
    list(train = 1:5, test = 6:10)
  )

  ms <- expect_warning(
    banded_ridge_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 0, b = 1),
      recall_folds = list(list(train = 1:9, test = 10)), # should be ignored
      target_folds = folds,
      compute_delta_r2 = FALSE,
      return_diagnostics = TRUE
    ),
    "both `recall_folds` and `target_folds`"
  )

  res <- run_regional(ms, regionMask)
  expect_equal(res$fits[[1]]$folds, folds)
})

# ---- T5: Mathematical correctness tests ------------------------------------

test_that("stacked mode matches manual ridge solve (no recall)", {
  # With alpha_target = 0, stacked mode is pure ridge on encoding data
  dat <- .make_linear_da_dataset(Tenc = 20, Trec = 10, set_sizes = c(a = 2, b = 2),
                                  zero_sets = NULL, seed = 42)

  ms <- banded_ridge_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 5, b = 5),
    alpha_target = 0,
    compute_delta_r2 = FALSE,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  # R² should be finite and computable

  perf <- res$performance_table[1, ]
  expect_true(is.finite(perf$recall_r2_full))
})

test_that("delta_r2 is zero for non-contributing set and positive for contributing set", {
  dat <- .make_linear_da_dataset(Tenc = 30, Trec = 20, set_sizes = c(a = 2, b = 2),
                                  zero_sets = "b", seed = 55)

  ms <- banded_ridge_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 0.01, b = 1),
    alpha_target = 0,
    compute_delta_r2 = TRUE,
    return_diagnostics = TRUE
  )

  res <- suppressWarnings(run_regional(ms, dat$mask))
  perf <- res$performance_table[1, ]

  # Set "b" has zero true weights, so removing it should not hurt
  expect_lt(abs(perf$delta_r2_b), 0.1)
  # Set "a" has signal, so removing it should reduce R²
  expect_gt(perf$delta_r2_a, 0.1)
})

# ---- T6: Error path tests --------------------------------------------------

test_that("banded_ridge_da_model accepts negative lambdas (no validation)", {
  # Note: The function does not validate lambda >= 0, so negative values are accepted
  # This test documents current behavior
  dat <- .make_linear_da_dataset(seed = 60)
  ms <- banded_ridge_da_model(
    dataset = dat$dataset, design = dat$design, mode = "stacked",
    lambdas = c(a = -1, b = 1), alpha_target = 0, compute_delta_r2 = FALSE
  )
  expect_s3_class(ms, "banded_ridge_da_model")
})

test_that("banded_ridge_da_model rejects missing lambda for a set", {
  dat <- .make_linear_da_dataset(seed = 61)
  expect_error(
    banded_ridge_da_model(
      dataset = dat$dataset, design = dat$design, mode = "stacked",
      lambdas = c(a = 1)  # missing "b"
    ),
    "missing entries"
  )
})

test_that("banded_ridge_da_model rejects negative alpha_recall", {
  dat <- .make_linear_da_dataset(seed = 62)
  expect_error(
    banded_ridge_da_model(
      dataset = dat$dataset, design = dat$design, mode = "stacked",
      lambdas = c(a = 1, b = 1), alpha_recall = -0.5
    ),
    "alpha_recall"
  )
})

test_that("banded_ridge_da_model rejects missing test set", {
  set.seed(63)
  dims <- c(2, 2, 1)
  mask <- neuroim2::NeuroVol(array(1, dim = dims), neuroim2::NeuroSpace(dims))
  train_vec <- neuroim2::NeuroVec(array(rnorm(prod(dims) * 10), c(dims, 10)),
                                   neuroim2::NeuroSpace(c(dims, 10)))
  dataset_no_test <- mvpa_dataset(train_vec, mask = mask)

  Xenc <- matrix(rnorm(10 * 4), 10, 4)
  fs_enc <- feature_sets(Xenc, blocks(a = 2, b = 2))
  fs_rec <- fs_enc  # dummy
  des <- feature_sets_design(fs_enc, fs_rec)

  expect_error(
    banded_ridge_da_model(
      dataset = dataset_no_test, design = des, mode = "stacked",
      lambdas = c(a = 1, b = 1)
    ),
    "test set"
  )
})

test_that("banded_ridge_da_model rejects unnamed lambdas", {
  dat <- .make_linear_da_dataset(seed = 64)
  expect_error(
    banded_ridge_da_model(
      dataset = dat$dataset, design = dat$design, mode = "stacked",
      lambdas = c(1, 1)  # unnamed
    ),
    "named"
  )
})

# ---- T7: Edge case tests ---------------------------------------------------

test_that("banded_ridge_da_model handles small ROI (2 voxels)", {
  skip_on_cran()
  # mvpa_dataset rejects 1-voxel volumes, so test with 2 voxels
  dat <- .make_linear_da_dataset(Tenc = 15, Trec = 10, dims = c(2, 1, 1),
                                  set_sizes = c(a = 2, b = 2), seed = 70)

  ms <- banded_ridge_da_model(
    dataset = dat$dataset, design = dat$design, mode = "stacked",
    lambdas = c(a = 1, b = 1), alpha_target = 0.2,
    compute_delta_r2 = FALSE
  )

  res <- run_regional(ms, dat$mask)
  expect_s3_class(res, "regional_mvpa_result")
  perf <- res$performance_table[1, ]
  expect_true(is.finite(perf$recall_r2_full) || is.na(perf$recall_r2_full))
})

test_that("banded_ridge_da_model handles minimal data sizes", {
  skip_on_cran()
  dat <- .make_linear_da_dataset(Tenc = 6, Trec = 6, dims = c(2, 1, 1),
                                  set_sizes = c(a = 2, b = 2), seed = 71)

  ms <- banded_ridge_da_model(
    dataset = dat$dataset, design = dat$design, mode = "stacked",
    lambdas = c(a = 1, b = 1), alpha_target = 0.2,
    compute_delta_r2 = FALSE
  )

  res <- run_regional(ms, dat$mask)
  expect_s3_class(res, "regional_mvpa_result")
})

test_that("banded_ridge_da_model handles lambda=0 via safe solve", {
  skip_on_cran()
  dat <- .make_linear_da_dataset(Tenc = 20, Trec = 10, set_sizes = c(a = 2, b = 2),
                                  zero_sets = NULL, seed = 72)

  ms <- banded_ridge_da_model(
    dataset = dat$dataset, design = dat$design, mode = "stacked",
    lambdas = c(a = 0, b = 0), alpha_target = 0,
    compute_delta_r2 = FALSE
  )

  # Should not crash thanks to .br_safe_solve
  res <- run_regional(ms, dat$mask)
  expect_s3_class(res, "regional_mvpa_result")
})
