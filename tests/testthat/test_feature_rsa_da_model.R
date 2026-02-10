context("feature_rsa_da_model")

.make_frda_dataset <- function(Tenc = 24, Trec = 18, dims = c(3, 3, 2), Dfeat = 8, seed = 101) {
  set.seed(seed)
  toy <- gen_sample_dataset(
    D = dims,
    nobs = Tenc,
    nlevels = 2,
    blocks = 2,
    external_test = TRUE,
    ntest_obs = Trec
  )

  Xenc <- matrix(rnorm(Tenc * Dfeat), Tenc, Dfeat)
  fs_enc <- feature_sets(Xenc, blocks(a = Dfeat / 2, b = Dfeat / 2))

  G <- matrix(runif(Trec * Tenc), Trec, Tenc)
  G <- G / rowSums(G)
  fs_rec <- expected_features(fs_enc, G, drop_null = FALSE, renormalize = TRUE)

  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, length.out = Trec))
  mask <- neuroim2::NeuroVol(array(1, dims), neuroim2::space(toy$dataset$mask))

  list(dataset = toy$dataset, design = des, mask = mask)
}

.make_frda_linear_dataset <- function(Tenc = 24, Trec = 18, dims = c(3, 3, 2), set_sizes = c(a = 4, b = 4), seed = 301) {
  set.seed(seed)
  Dfeat <- sum(set_sizes)
  V <- prod(dims)

  X_enc <- matrix(rnorm(Tenc * Dfeat), Tenc, Dfeat)
  X_rec <- matrix(rnorm(Trec * Dfeat), Trec, Dfeat)
  B <- matrix(rnorm(Dfeat * V), nrow = Dfeat, ncol = V)

  Y_enc <- X_enc %*% B
  Y_rec <- X_rec %*% B

  arr_enc <- array(0, dim = c(dims, Tenc))
  for (t in seq_len(Tenc)) arr_enc[,,, t] <- array(Y_enc[t, ], dim = dims)
  arr_rec <- array(0, dim = c(dims, Trec))
  for (t in seq_len(Trec)) arr_rec[,,, t] <- array(Y_rec[t, ], dim = dims)

  mask <- neuroim2::NeuroVol(array(1, dim = dims), neuroim2::NeuroSpace(dims))
  train_vec <- neuroim2::NeuroVec(arr_enc, neuroim2::NeuroSpace(c(dims, Tenc)))
  test_vec <- neuroim2::NeuroVec(arr_rec, neuroim2::NeuroSpace(c(dims, Trec)))
  dataset <- mvpa_dataset(train_vec, test_vec, mask)

  spec <- do.call(blocks, as.list(set_sizes))
  fs_enc <- feature_sets(X_enc, spec)
  fs_rec <- feature_sets(X_rec, spec, set_order = fs_enc$set_order)
  des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, length.out = Trec))

  list(dataset = dataset, design = des, mask = mask)
}

test_that("feature_rsa_da_model runs and returns target geometry metrics", {
  skip_on_cran()
  dat <- .make_frda_dataset(seed = 201)

  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  expect_s3_class(res, "regional_mvpa_result")
  perf_names <- names(res$performance_table)

  expect_true("target_pattern_correlation" %in% perf_names)
  expect_true("target_pattern_discrimination" %in% perf_names)
  expect_true("target_rdm_correlation" %in% perf_names)
  expect_true("target_voxel_correlation" %in% perf_names)
  expect_true("target_mse_full" %in% perf_names)
  expect_true("target_r2_full" %in% perf_names)
})

test_that("feature_rsa_da_model recovers exact mapping in stacked mode", {
  dat <- .make_frda_linear_dataset(seed = 302, Tenc = 28, Trec = 20)
  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 0, b = 0),
    alpha_target = 0,
    rsa_simfun = "spearman"
  )

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_gt(perf$target_r2_full, 0.99)
  expect_gt(perf$target_rdm_correlation, 0.95)
  expect_gt(perf$target_pattern_correlation, 0.90)
})

test_that("feature_rsa_da_model validates lambda names", {
  dat <- .make_frda_dataset(seed = 303)
  expect_error(
    feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1)
    ),
    "missing entries"
  )
})

test_that("feature_rsa_da_model validates non-negative alpha_target", {
  dat <- .make_frda_dataset(seed = 304)
  expect_error(
    feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1, b = 1),
      alpha_target = -0.2
    ),
    "alpha_target"
  )
})

test_that("feature_rsa_da_model supports single-run purge gap", {
  skip_on_cran()
  dat <- .make_frda_dataset(seed = 202, Trec = 18)
  dat$design$block_var_test <- rep(1, 18)

  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    recall_nfolds = 3,
    target_gap = 1,
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  folds <- res$fits[[1]]$folds
  expect_equal(length(folds), 3)

  for (f in folds) {
    test_idx <- f$test
    purge_lo <- max(1L, min(test_idx) - 1L)
    purge_hi <- min(18L, max(test_idx) + 1L)
    purged <- seq.int(purge_lo, purge_hi)
    expect_false(any(f$train %in% purged))
  }
})

test_that("feature_rsa_da_model computes single-run permutation p-values", {
  skip_on_cran()
  dat <- .make_frda_dataset(seed = 203, Trec = 16)
  dat$design$block_var_test <- rep(1, 16)

  ms <- feature_rsa_da_model(
    dataset = dat$dataset,
    design = dat$design,
    mode = "stacked",
    lambdas = c(a = 1, b = 1),
    recall_nfolds = 4,
    target_nperm = 6,
    target_perm_strategy = "circular_shift",
    return_diagnostics = TRUE
  )

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_true("target_perm_p_rdm_correlation" %in% names(perf))
  expect_true("target_perm_p_mse_full" %in% names(perf))
  expect_true(perf$target_perm_p_rdm_correlation >= 0 && perf$target_perm_p_rdm_correlation <= 1)
  expect_true(perf$target_perm_p_mse_full >= 0 && perf$target_perm_p_mse_full <= 1)
  expect_equal(perf$target_perm_n, 6)
})

test_that("feature_rsa_da_model disables permutation null for multi-run targets", {
  dat <- .make_frda_dataset(seed = 204, Trec = 18)
  dat$design$block_var_test <- rep(1:2, each = 9)

  ms <- expect_warning(
    feature_rsa_da_model(
      dataset = dat$dataset,
      design = dat$design,
      mode = "stacked",
      lambdas = c(a = 1, b = 1),
      target_nperm = 5
    ),
    "single-run targets"
  )
  expect_equal(ms$target_nperm, 0L)

  res <- run_regional(ms, dat$mask)
  perf <- res$performance_table[1, ]
  expect_false("target_perm_p_rdm_correlation" %in% names(perf))
})
