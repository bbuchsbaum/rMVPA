context("era_partition_model")

.era_partition_test_model <- function(K = 10, p = 4, include_procrustes = TRUE,
                                      first_order_nuisance = NULL,
                                      second_order_nuisance = NULL,
                                      item_block_enc = NULL,
                                      item_block_ret = NULL,
                                      item_run_enc = NULL,
                                      item_run_ret = NULL,
                                      item_time_enc = NULL,
                                      item_time_ret = NULL,
                                      item_category = NULL,
                                      auto_nuisance = TRUE,
                                      global_nuisance = FALSE) {
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
    item_run_enc = item_run_enc,
    item_run_ret = item_run_ret,
    item_time_enc = item_time_enc,
    item_time_ret = item_time_ret,
    item_category = item_category,
    include_procrustes = include_procrustes,
    min_procrustes_train_items = 3L,
    auto_nuisance = auto_nuisance,
    global_nuisance = global_nuisance
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

.era_partition_xdec_fixture <- function(K = 12, p = 6, seed = 100) {
  set.seed(seed)
  toy <- gen_sample_dataset(
    D = c(3, 3, 3),
    nobs = K,
    nlevels = 3,
    blocks = 3,
    external_test = TRUE,
    ntest_obs = K
  )
  keys <- paste0("item", seq_len(K))
  toy$design$train_design$item <- factor(keys, levels = keys)
  toy$design$test_design$item <- factor(keys, levels = keys)
  toy$design$train_design$category <- toy$design$y_train
  toy$design$test_design$category <- toy$design$y_test

  list(
    toy = toy,
    E = matrix(rnorm(K * p), K, p),
    R = matrix(rnorm(K * p), K, p),
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
  run_enc <- setNames(c("r1", "r2", "r1", "r2"), keys)
  run_ret <- setNames(c("r1", "r1", "r2", "r2"), keys)
  time_enc <- setNames(c(10, 20, 30, 40), keys)
  time_ret <- setNames(c(1, 2, 3, 4), keys)

  model <- .era_partition_test_model(
    K = K,
    p = 4,
    item_block_enc = block_enc,
    item_block_ret = block_ret,
    item_run_enc = run_enc,
    item_run_ret = run_ret,
    item_time_enc = time_enc,
    item_time_ret = time_ret,
    include_procrustes = FALSE
  )

  first <- rMVPA:::.era_partition_first_nuisance(model, keys)
  expect_equal(unname(matrix(first$same_block_cross, K, K)), unname(outer(block_ret, block_enc, "==") * 1))
  expect_equal(unname(matrix(first$same_run_cross, K, K)), unname(outer(run_ret, run_enc, "==") * 1))
  expect_equal(matrix(first$enc_time, K, K), matrix(rep(time_enc, each = K), K, K))
  expect_equal(matrix(first$ret_time, K, K), matrix(rep(time_ret, times = K), K, K))
  expect_equal(unname(matrix(first$abs_lag, K, K)), unname(abs(outer(time_ret, time_enc, "-"))))

  second <- rMVPA:::.era_partition_second_nuisance(model, keys)
  run_enc_mat <- outer(run_enc, run_enc, "==")
  run_ret_mat <- outer(run_ret, run_ret, "==")
  expect_equal(second$same_run_enc, as.numeric(run_enc_mat[lower.tri(run_enc_mat)]))
  expect_equal(second$same_run_ret, as.numeric(run_ret_mat[lower.tri(run_ret_mat)]))
  enc_td <- abs(outer(time_enc, time_enc, "-"))
  ret_td <- abs(outer(time_ret, time_ret, "-"))
  expect_equal(second$temporal_distance_enc, enc_td[lower.tri(enc_td)])
  expect_equal(second$temporal_distance_ret, ret_td[lower.tri(ret_td)])
})

test_that("era_partition_model can disable or select automatic nuisance groups", {
  K <- 5
  keys <- paste0("item", seq_len(K))
  block_enc <- setNames(rep(c("b1", "b2"), length.out = K), keys)
  block_ret <- setNames(rep(c("b2", "b1"), length.out = K), keys)
  run_enc <- setNames(rep(c("r1", "r2"), length.out = K), keys)
  run_ret <- setNames(rep(c("s1", "s2"), length.out = K), keys)
  time_enc <- setNames(seq_len(K), keys)
  time_ret <- setNames(seq_len(K) + 10, keys)
  category <- setNames(rep(letters[1:2], length.out = K), keys)

  model_none <- .era_partition_test_model(
    K = K, p = 3,
    item_block_enc = block_enc,
    item_block_ret = block_ret,
    item_run_enc = run_enc,
    item_run_ret = run_ret,
    item_time_enc = time_enc,
    item_time_ret = time_ret,
    item_category = category,
    include_procrustes = FALSE,
    auto_nuisance = FALSE
  )
  expect_length(rMVPA:::.era_partition_first_nuisance(model_none, keys), 0L)
  expect_length(rMVPA:::.era_partition_second_nuisance(model_none, keys), 0L)

  model_run <- .era_partition_test_model(
    K = K, p = 3,
    item_block_enc = block_enc,
    item_block_ret = block_ret,
    item_run_enc = run_enc,
    item_run_ret = run_ret,
    item_time_enc = time_enc,
    item_time_ret = time_ret,
    item_category = category,
    include_procrustes = FALSE,
    auto_nuisance = "run"
  )
  expect_identical(names(rMVPA:::.era_partition_first_nuisance(model_run, keys)), "same_run_cross")
  expect_identical(names(rMVPA:::.era_partition_second_nuisance(model_run, keys)), c("same_run_enc", "same_run_ret"))
})

test_that("era_partition_model consumes global nuisance as a selectable auto group", {
  K <- 4
  keys <- paste0("item", seq_len(K))
  S_cross <- matrix(seq_len(K * K), K, K)
  D_enc <- as.matrix(stats::dist(matrix(c(0, 1, 3, 6), ncol = 1)))
  D_ret <- as.matrix(stats::dist(matrix(c(0, 2, 5, 9), ncol = 1)))
  dimnames(S_cross) <- list(keys, keys)
  rownames(D_enc) <- colnames(D_enc) <- keys
  rownames(D_ret) <- colnames(D_ret) <- keys
  global <- list(S_cross = S_cross, D_enc = D_enc, D_ret = D_ret)

  model_global <- .era_partition_test_model(
    K = K, p = 3,
    include_procrustes = FALSE,
    auto_nuisance = "global",
    global_nuisance = global
  )
  first <- rMVPA:::.era_partition_first_nuisance(model_global, keys)
  second <- rMVPA:::.era_partition_second_nuisance(model_global, keys)

  expect_identical(names(first), "global_cross")
  expect_equal(first$global_cross, as.numeric(S_cross))
  expect_identical(names(second), c("global_enc", "global_ret"))
  expect_equal(second$global_enc, as.numeric(D_enc[lower.tri(D_enc)]))
  expect_equal(second$global_ret, as.numeric(D_ret[lower.tri(D_ret)]))

  model_off <- .era_partition_test_model(
    K = K, p = 3,
    include_procrustes = FALSE,
    auto_nuisance = FALSE,
    global_nuisance = global
  )
  expect_length(rMVPA:::.era_partition_first_nuisance(model_off, keys), 0L)
  expect_length(rMVPA:::.era_partition_second_nuisance(model_off, keys), 0L)
})

test_that("era_partition_model compares phase-scoped factor run labels", {
  K <- 4
  keys <- paste0("item", seq_len(K))
  run_enc <- factor(setNames(c("enc_1", "enc_2", "enc_1", "enc_2"), keys))
  run_ret <- factor(setNames(c("ret_1", "ret_1", "ret_2", "ret_2"), keys))

  model <- .era_partition_test_model(
    K = K, p = 3,
    item_run_enc = run_enc,
    item_run_ret = run_ret,
    include_procrustes = FALSE,
    auto_nuisance = "run"
  )

  first <- rMVPA:::.era_partition_first_nuisance(model, keys)
  second <- rMVPA:::.era_partition_second_nuisance(model, keys)

  expect_equal(unname(matrix(first$same_run_cross, K, K)), matrix(0, K, K))
  expect_equal(second$same_run_enc, as.numeric(outer(as.character(run_enc), as.character(run_enc), "==")[lower.tri(diag(K))]))
  expect_equal(second$same_run_ret, as.numeric(outer(as.character(run_ret), as.character(run_ret), "==")[lower.tri(diag(K))]))
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

test_that("era_partition_model direct xdec metrics match naive_xdec_model", {
  fix <- .era_partition_xdec_fixture(K = 36, p = 27, seed = 10)
  toy <- fix$toy

  roi_data <- list(
    train_data = t(as.matrix(toy$dataset$train_data)),
    test_data = t(as.matrix(toy$dataset$test_data)),
    indices = seq_len(nrow(as.matrix(toy$dataset$train_data)))
  )
  context <- fix$context

  era_model <- era_partition_model(
    dataset = toy$dataset,
    design = toy$design,
    key_var = ~ item,
    distfun = eucdist(),
    xdec_link_by = "category",
    include_procrustes = FALSE
  )
  naive_model <- naive_xdec_model(
    dataset = toy$dataset,
    design = toy$design,
    link_by = "category",
    return_predictions = FALSE
  )

  era_out <- fit_roi(era_model, roi_data, context)
  naive_out <- fit_roi(naive_model, roi_data, context)

  expect_false(era_out$error)
  expect_false(naive_out$error)
  expect_equal(
    era_out$metrics[c("xdec_Accuracy", "xdec_AUC")],
    setNames(naive_out$metrics[c("Accuracy", "AUC")], c("xdec_Accuracy", "xdec_AUC")),
    tolerance = 1e-12
  )
})

test_that("era_partition_model xdec schema follows compute flag", {
  fix <- .era_partition_xdec_fixture()

  with_xdec <- era_partition_model(
    dataset = fix$toy$dataset,
    design = fix$toy$design,
    key_var = ~ item,
    xdec_link_by = "category",
    include_procrustes = FALSE
  )
  without_xdec <- era_partition_model(
    dataset = fix$toy$dataset,
    design = fix$toy$design,
    key_var = ~ item,
    compute_xdec_performance = FALSE,
    include_procrustes = FALSE
  )

  expect_true(all(c("xdec_Accuracy", "xdec_AUC") %in% names(output_schema(with_xdec))))
  expect_false(any(grepl("^xdec_", names(output_schema(without_xdec)))))

  out <- .era_partition_fit_direct(without_xdec, fix$E, fix$R)
  expect_false(out$error)
  expect_false(any(grepl("^xdec_", names(out$metrics))))
})

test_that("era_partition_model stores diagnostic matrices and xdec predictions on request", {
  fix <- .era_partition_xdec_fixture(K = 15, p = 8, seed = 101)
  model <- era_partition_model(
    dataset = fix$toy$dataset,
    design = fix$toy$design,
    key_var = ~ item,
    distfun = eucdist(),
    xdec_link_by = "category",
    include_procrustes = FALSE,
    return_matrices = TRUE,
    return_xdec_predictions = TRUE
  )

  out <- .era_partition_fit_direct(model, fix$E, fix$R)

  expect_false(out$error)
  expect_true(all(c(
    "keys",
    "source_prototypes",
    "target_prototypes",
    "cross_similarity",
    "source_rdm",
    "target_rdm",
    "procrustes_similarity",
    "xdec_result"
  ) %in% names(out$result)))
  expect_equal(length(out$result$keys), 15)
  expect_equal(dim(out$result$cross_similarity), c(15L, 15L))
  expect_equal(rownames(out$result$cross_similarity), out$result$keys)
  expect_equal(colnames(out$result$cross_similarity), out$result$keys)
  expect_true(all(is.finite(diag(out$result$source_rdm))))
  expect_equal(unname(diag(out$result$source_rdm)), rep(0, 15), tolerance = 1e-12)
  expect_equal(unname(diag(out$result$target_rdm)), rep(0, 15), tolerance = 1e-12)

  probs <- out$result$xdec_result$probs
  expect_true(all(is.finite(unname(probs))))
  expect_equal(rowSums(probs), rep(1, nrow(probs)), tolerance = 1e-12)
  expect_equal(nrow(probs), nrow(fix$R))
  expect_true(all(c("xdec_Accuracy", "xdec_AUC") %in% names(out$metrics)))
})

test_that("era_partition_model falls back to key_var when requested xdec link is incomplete", {
  fix <- .era_partition_xdec_fixture(K = 10, p = 5, seed = 102)
  fix$toy$design$train_design$train_only <- fix$toy$design$train_design$category

  model <- era_partition_model(
    dataset = fix$toy$dataset,
    design = fix$toy$design,
    key_var = ~ item,
    xdec_link_by = "train_only",
    include_procrustes = FALSE,
    return_xdec_predictions = TRUE
  )

  out <- .era_partition_fit_direct(model, fix$E, fix$R)

  expect_null(model$xdec_link_by)
  expect_equal(model$xdec_link_by_requested, "train_only")
  expect_false(out$error)
  expect_equal(ncol(out$result$xdec_result$probs), 10)
  expect_equal(colnames(out$result$xdec_result$probs), levels(fix$toy$design$train_design$item))
  expect_true(all(is.finite(out$metrics[c("xdec_Accuracy", "xdec_AUC")])))
})

test_that("era_partition_model reports NA xdec metrics for degenerate xdec labels only", {
  fix <- .era_partition_xdec_fixture(K = 8, p = 5, seed = 103)
  fix$toy$design$train_design$single_class <- "same"
  fix$toy$design$test_design$single_class <- "same"

  model <- era_partition_model(
    dataset = fix$toy$dataset,
    design = fix$toy$design,
    key_var = ~ item,
    xdec_link_by = "single_class",
    include_procrustes = FALSE
  )

  out <- .era_partition_fit_direct(model, fix$E, fix$R)

  expect_false(out$error)
  expect_true(all(c("xdec_Accuracy", "xdec_AUC") %in% names(out$metrics)))
  expect_true(all(is.na(out$metrics[c("xdec_Accuracy", "xdec_AUC")])))
  expect_true(is.finite(out$metrics[["naive_top1_acc"]]))
  expect_true(is.finite(out$metrics[["geom_cor"]]))
})

test_that("era_partition_model direct xdec is invariant to paired row order", {
  fix <- .era_partition_xdec_fixture(K = 18, p = 7, seed = 104)
  model <- era_partition_model(
    dataset = fix$toy$dataset,
    design = fix$toy$design,
    key_var = ~ item,
    xdec_link_by = "category",
    include_procrustes = FALSE
  )
  ref <- .era_partition_fit_direct(model, fix$E, fix$R)

  perm <- sample(seq_len(nrow(fix$E)))
  model_perm <- model
  model_perm$design$train_design <- model$design$train_design[perm, , drop = FALSE]
  model_perm$design$test_design <- model$design$test_design[perm, , drop = FALSE]
  model_perm$design$y_train <- model$design$y_train[perm]
  model_perm$design$y_test <- model$design$y_test[perm]
  permuted <- .era_partition_fit_direct(
    model_perm,
    fix$E[perm, , drop = FALSE],
    fix$R[perm, , drop = FALSE]
  )

  invariant <- c(
    "xdec_Accuracy",
    "xdec_AUC",
    "naive_top1_acc",
    "first_order_delta_r2",
    "second_order_delta_r2",
    "geom_cor"
  )
  expect_false(ref$error)
  expect_false(permuted$error)
  expect_equal(permuted$metrics[invariant], ref$metrics[invariant], tolerance = 1e-10)
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
  run_enc <- setNames(rep(c("enc_1", "enc_2"), length.out = K), keys)
  run_ret <- setNames(rep(c("ret_1", "ret_2"), length.out = K), keys)

  model <- .era_partition_test_model(
    K,
    p,
    item_block_enc = block_enc,
    item_block_ret = block_ret,
    item_run_enc = run_enc,
    item_run_ret = run_ret,
    item_time_enc = time_enc,
    item_time_ret = time_ret,
    item_category = category
  )
  out <- .era_partition_fit_direct(model, E, R)

  expect_false(out$error)
  expect_equal(out$metrics[["nuisance_first_order_n"]], 6)
  expect_equal(out$metrics[["nuisance_second_order_n"]], 7)
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

test_that("era_partition_model representative ROI fit stays within perf budget", {
  skip_on_cran()
  skip_if_not_perf_tests()

  fix <- .era_partition_xdec_fixture(K = 60, p = 48, seed = 105)
  model <- era_partition_model(
    dataset = fix$toy$dataset,
    design = fix$toy$design,
    key_var = ~ item,
    distfun = eucdist(),
    xdec_link_by = "category",
    include_procrustes = TRUE,
    min_procrustes_train_items = 5L
  )

  elapsed <- as.numeric(system.time({
    out <- .era_partition_fit_direct(model, fix$E, fix$R)
  })["elapsed"])
  budget <- as.numeric(Sys.getenv("RMVPA_ERA_PARTITION_MAX_SECONDS", "2.5"))

  message(sprintf(
    "era_partition perf guardrail: elapsed=%.3fs budget=%.3fs K=%d p=%d",
    elapsed, budget, nrow(fix$E), ncol(fix$E)
  ))
  expect_false(out$error)
  expect_lte(elapsed, budget)
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

test_that("era_partition_model warns when item-level block metadata is missing", {
  expect_warning(
    .era_partition_test_model(K = 8, p = 3, include_procrustes = FALSE),
    "item_block_enc.*item_block_ret"
  )

  toy <- gen_sample_dataset(D = c(3,3,3), nobs = 8, nlevels = 2, blocks = 2,
                            external_test = TRUE, ntest_obs = 8)
  keys <- paste0("item", seq_len(8))
  toy$design$train_design$item <- factor(keys, levels = keys)
  toy$design$test_design$item  <- factor(keys, levels = keys)
  expect_error(
    era_partition_model(
      dataset = toy$dataset, design = toy$design,
      key_var = ~ item, distfun = eucdist(),
      include_procrustes = FALSE,
      require_run_metadata = TRUE
    ),
    "item_block_enc"
  )
})

test_that("era_partition_model flags overlapping enc/ret block labels", {
  K <- 6
  keys <- paste0("item", seq_len(K))
  block_enc <- setNames(rep(c(1, 2), length.out = K), keys)
  block_ret <- setNames(rep(c(1, 2), length.out = K), keys)
  expect_warning(
    .era_partition_test_model(
      K = K, p = 3,
      item_block_enc = block_enc,
      item_block_ret = block_ret,
      include_procrustes = FALSE
    ),
    "phase-scoped labels"
  )
})
