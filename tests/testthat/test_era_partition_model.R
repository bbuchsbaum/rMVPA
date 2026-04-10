context("era_partition_model")

.era_partition_test_model <- function(K = 10, p = 4, include_procrustes = TRUE,
                                      first_order_nuisance = NULL,
                                      second_order_nuisance = NULL,
                                      item_block_enc = NULL,
                                      item_block_ret = NULL,
                                      item_time_enc = NULL,
                                      item_time_ret = NULL,
                                      item_category = NULL) {
  toy <- gen_sample_dataset(
    D = c(3, 3, 3),
    nobs = K,
    nlevels = 2,
    blocks = 2,
    external_test = TRUE,
    ntest_obs = K
  )
  keys <- paste0("item", seq_len(K))
  toy$design$train_design$item <- factor(keys, levels = keys)
  toy$design$test_design$item <- factor(keys, levels = keys)

  era_partition_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = ~ item,
    distfun = eucdist(),
    first_order_nuisance = first_order_nuisance,
    second_order_nuisance = second_order_nuisance,
    item_block_enc = item_block_enc,
    item_block_ret = item_block_ret,
    item_time_enc = item_time_enc,
    item_time_ret = item_time_ret,
    item_category = item_category,
    include_procrustes = include_procrustes,
    min_procrustes_train_items = 3L
  )
}

.era_partition_fit_direct <- function(model, E, R) {
  fit_roi(
    model,
    roi_data = list(
      train_data = E,
      test_data = R,
      indices = seq_len(ncol(E))
    ),
    context = list(id = 1L)
  )
}

test_that("delta R2 matches an explicit nested-lm oracle", {
  y <- c(1, 2, 1.5, 4, 3.5, 3, 5, 5.5)
  signal <- c(0, 0.2, 0.3, 1, 1.1, 0.9, 1.8, 2)
  nuisance <- list(
    drift = seq_along(y),
    block = rep(c(0, 1), length.out = length(y))
  )

  got <- rMVPA:::.era_partition_delta_r2(y, signal, nuisance)

  dat <- data.frame(y = y, signal = signal, drift = nuisance$drift, block = nuisance$block)
  full <- lm(y ~ signal + drift + block, data = dat)
  nuis <- lm(y ~ drift + block, data = dat)
  r2_full <- 1 - sum(resid(full)^2) / sum((y - mean(y))^2)
  r2_nuis <- 1 - sum(resid(nuis)^2) / sum((y - mean(y))^2)

  expect_equal(got$full_r2, r2_full, tolerance = 1e-12)
  expect_equal(got$nuisance_r2, r2_nuis, tolerance = 1e-12)
  expect_equal(got$delta_r2, r2_full - r2_nuis, tolerance = 1e-12)
  expect_equal(got$partial_r2, (r2_full - r2_nuis) / (1 - r2_nuis), tolerance = 1e-12)
})

test_that("nuisance builders preserve first-order and second-order orientation", {
  K <- 4
  keys <- paste0("item", seq_len(K))
  block_enc <- setNames(c("a", "b", "a", "b"), keys)
  block_ret <- setNames(c("a", "a", "b", "b"), keys)
  time_enc <- setNames(c(10, 20, 30, 40), keys)
  time_ret <- setNames(c(1, 2, 3, 4), keys)

  model <- .era_partition_test_model(
    K = K,
    p = 4,
    item_block_enc = block_enc,
    item_block_ret = block_ret,
    item_time_enc = time_enc,
    item_time_ret = time_ret,
    include_procrustes = FALSE
  )

  first <- rMVPA:::.era_partition_first_nuisance(model, keys)
  expect_equal(unname(matrix(first$same_block_cross, K, K)), unname(outer(block_ret, block_enc, "==") * 1))
  expect_equal(matrix(first$enc_time, K, K), matrix(rep(time_enc, each = K), K, K))
  expect_equal(matrix(first$ret_time, K, K), matrix(rep(time_ret, times = K), K, K))
  expect_equal(unname(matrix(first$abs_lag, K, K)), unname(abs(outer(time_ret, time_enc, "-"))))

  second <- rMVPA:::.era_partition_second_nuisance(model, keys)
  enc_td <- abs(outer(time_enc, time_enc, "-"))
  ret_td <- abs(outer(time_ret, time_ret, "-"))
  expect_equal(second$temporal_distance_enc, enc_td[lower.tri(enc_td)])
  expect_equal(second$temporal_distance_ret, ret_td[lower.tri(ret_td)])
})

test_that("era_partition_model returns matched first- and second-order metrics", {
  set.seed(1)

  K <- 8
  p <- 30
  E <- matrix(rnorm(K * p), K, p)
  R <- E + matrix(rnorm(K * p, sd = 0.001), K, p)

  model <- .era_partition_test_model(K, p, include_procrustes = FALSE)
  out <- .era_partition_fit_direct(model, E, R)

  expect_false(out$error)
  expect_true(all(c(
    "first_order_delta_r2",
    "second_order_delta_r2",
    "naive_top1_acc",
    "geom_cor",
    "procrustes_top1_acc"
  ) %in% names(out$metrics)))
  expect_gt(out$metrics[["first_order_delta_r2"]], 0.75)
  expect_gt(out$metrics[["second_order_delta_r2"]], 0.9)
  expect_equal(out$metrics[["naive_top1_acc"]], 1)
})

test_that("era_partition_model is invariant to paired item row permutations", {
  set.seed(11)

  K <- 9
  p <- 12
  keys <- paste0("item", seq_len(K))
  E <- matrix(rnorm(K * p), K, p)
  R <- E + matrix(rnorm(K * p, sd = 0.01), K, p)

  model_a <- .era_partition_test_model(K, p, include_procrustes = FALSE)
  out_a <- .era_partition_fit_direct(model_a, E, R)

  perm <- sample(seq_len(K))
  model_b <- .era_partition_test_model(K, p, include_procrustes = FALSE)
  model_b$design$train_design$item <- factor(keys[perm], levels = keys)
  model_b$design$test_design$item <- factor(keys[perm], levels = keys)
  out_b <- .era_partition_fit_direct(model_b, E[perm, , drop = FALSE], R[perm, , drop = FALSE])

  invariant <- c(
    "first_order_delta_r2",
    "second_order_delta_r2",
    "naive_top1_acc",
    "geom_cor"
  )
  expect_equal(out_b$metrics[invariant], out_a$metrics[invariant], tolerance = 1e-10)
})

test_that("era_partition_model detects preserved geometry under orthogonal rotation", {
  set.seed(2)

  K <- 12
  p <- 4
  E <- matrix(rnorm(K * p), K, p)
  q <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  R <- E %*% q + matrix(rnorm(K * p, sd = 0.005), K, p)

  model <- .era_partition_test_model(K, p)
  out <- .era_partition_fit_direct(model, E, R)

  expect_false(out$error)
  expect_gt(out$metrics[["second_order_delta_r2"]], 0.95)
  expect_gt(out$metrics[["geom_cor"]], 0.95)
  expect_gt(out$metrics[["procrustes_top1_acc"]], 0.8)
})

test_that("era_partition_model includes automatic block and temporal nuisances", {
  set.seed(3)

  K <- 10
  p <- 4
  keys <- paste0("item", seq_len(K))
  E <- matrix(rnorm(K * p), K, p)
  R <- E + matrix(rnorm(K * p, sd = 0.02), K, p)

  block_enc <- setNames(rep(1:2, length.out = K), keys)
  block_ret <- setNames(rep(2:1, length.out = K), keys)
  time_enc <- setNames(seq_len(K), keys)
  time_ret <- setNames(seq_len(K) + 10, keys)
  category <- setNames(rep(letters[1:2], length.out = K), keys)

  model <- .era_partition_test_model(
    K,
    p,
    item_block_enc = block_enc,
    item_block_ret = block_ret,
    item_time_enc = time_enc,
    item_time_ret = time_ret,
    item_category = category
  )
  out <- .era_partition_fit_direct(model, E, R)

  expect_false(out$error)
  expect_equal(out$metrics[["nuisance_first_order_n"]], 5)
  expect_equal(out$metrics[["nuisance_second_order_n"]], 5)
  expect_true(is.finite(out$metrics[["first_order_delta_r2"]]))
  expect_true(is.finite(out$metrics[["second_order_delta_r2"]]))
})

test_that("Procrustes fit excludes the held-out target item", {
  set.seed(4)

  K <- 8
  p <- 4
  E <- matrix(rnorm(K * p), K, p)
  q <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
  R <- E %*% q

  holdout <- 3L
  train_idx <- setdiff(seq_len(K), holdout)
  fit_a <- rMVPA:::.era_partition_procrustes_fit(E, R, train_idx, center = TRUE)

  R_changed <- R
  R_changed[holdout, ] <- R_changed[holdout, ] + 100
  fit_b <- rMVPA:::.era_partition_procrustes_fit(E, R_changed, train_idx, center = TRUE)

  expect_equal(fit_a$T, fit_b$T, tolerance = 1e-10)
  expect_equal(fit_a$mu_R, fit_b$mu_R, tolerance = 1e-10)
})

test_that("degenerate variance-partition inputs return NA rather than misleading finite effects", {
  y <- rep(1, 8)
  signal <- rep(c(0, 1), 4)
  got_constant_y <- rMVPA:::.era_partition_delta_r2(y, signal)

  expect_true(is.na(got_constant_y$delta_r2))
  expect_equal(got_constant_y$n_obs, length(y))

  got_constant_signal <- rMVPA:::.era_partition_delta_r2(rnorm(8), rep(1, 8))
  expect_true(is.na(got_constant_signal$delta_r2))

  got_bad_nuisance <- try(
    rMVPA:::.era_partition_delta_r2(rnorm(8), signal, list(short = 1:3)),
    silent = TRUE
  )
  expect_s3_class(got_bad_nuisance, "try-error")
})

test_that("Procrustes scores return NA when leakage-free training set is too small", {
  set.seed(12)

  E <- matrix(rnorm(5 * 4), 5, 4)
  R <- matrix(rnorm(5 * 4), 5, 4)
  got <- rMVPA:::.era_partition_procrustes_scores(E, R, min_train_items = 5L)

  expect_true(all(is.na(got$scores)))
  expect_equal(got$train_n, rep(4, 5))
})

test_that("era_partition_model can run through regional iterator", {
  skip_on_cran()
  set.seed(5)

  toy <- gen_sample_dataset(D = c(4, 4, 4), nobs = 40, nlevels = 4,
                            blocks = 4, external_test = TRUE)
  toy$design$train_design$item <- toy$design$train_design$Y
  toy$design$test_design$item <- toy$design$test_design$Ytest

  region_mask <- neuroim2::NeuroVol(
    sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
    neuroim2::space(toy$dataset$mask)
  )

  model <- era_partition_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = ~ item
  )

  res <- run_regional(model, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(all(c("first_order_delta_r2", "second_order_delta_r2") %in% names(res$performance_table)))
})
