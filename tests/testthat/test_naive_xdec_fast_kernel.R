build_naive_xdec_fast_fixture <- function(
    nobs = 36L,
    ntest_obs = 36L,
    nlevels = 3L,
    blocks = 6L,
    p = 48L,
    seed = 9401L
) {
  set.seed(seed)
  ds <- gen_sample_dataset(
    D = c(5, 5, 5),
    nobs = nobs,
    nlevels = nlevels,
    blocks = blocks,
    external_test = TRUE,
    ntest_obs = ntest_obs
  )

  roi_data <- list(
    train_data = matrix(rnorm(nobs * p), nrow = nobs, ncol = p),
    test_data = matrix(rnorm(ntest_obs * p), nrow = ntest_obs, ncol = p),
    indices = seq_len(p)
  )

  list(dataset = ds$dataset, design = ds$design, roi_data = roi_data, context = list(id = 1L))
}

test_that("naive_xdec fast kernel cache is constructed by default", {
  fix <- build_naive_xdec_fast_fixture()

  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = FALSE)
  expect_true(is.list(model$.fast_kernel))
  expect_true(is.matrix(model$.fast_kernel$group_mat))
  expect_true(is.numeric(model$.fast_kernel$counts))
})

test_that("naive_xdec fast kernel preserves fit_roi outputs (differential parity)", {
  fix <- build_naive_xdec_fast_fixture(seed = 9402L)

  old_opt <- options(
    rMVPA.searchlight_mode = "legacy",
    rMVPA.naive_xdec_fast_kernel = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  base_model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = TRUE)
  base <- fit_roi(base_model, fix$roi_data, fix$context)
  expect_false(base$error)

  options(rMVPA.searchlight_mode = "fast", rMVPA.naive_xdec_fast_kernel = NULL)
  fast_model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = TRUE)
  fast <- fit_roi(fast_model, fix$roi_data, fix$context)
  expect_false(fast$error)

  expect_equal(base$metrics, fast$metrics, tolerance = 1e-10)
  expect_identical(as.character(base$result$predicted), as.character(fast$result$predicted))
  expect_equal(unname(base$result$probs), unname(fast$result$probs), tolerance = 1e-10)
})

test_that("naive_xdec fast kernel is invariant to feature permutation (metamorphic)", {
  fix <- build_naive_xdec_fast_fixture(seed = 9403L)
  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.naive_xdec_fast_kernel = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = TRUE)
  ref <- fit_roi(model, fix$roi_data, fix$context)
  expect_false(ref$error)

  perm <- sample.int(ncol(fix$roi_data$train_data))
  roi_perm <- list(
    train_data = fix$roi_data$train_data[, perm, drop = FALSE],
    test_data = fix$roi_data$test_data[, perm, drop = FALSE],
    indices = fix$roi_data$indices
  )
  permuted <- fit_roi(model, roi_perm, fix$context)
  expect_false(permuted$error)

  expect_equal(ref$metrics, permuted$metrics, tolerance = 1e-10)
  expect_identical(as.character(ref$result$predicted), as.character(permuted$result$predicted))
  expect_equal(unname(ref$result$probs), unname(permuted$result$probs), tolerance = 1e-10)
})

test_that("naive_xdec fast kernel handles class imbalance and unseen target labels", {
  fix <- build_naive_xdec_fast_fixture(seed = 9404L)

  y_tr <- as.character(fix$design$y_train)
  y_tr[] <- "a"
  y_tr[seq(1L, length(y_tr), by = 4L)] <- "b"
  y_tr[seq(2L, length(y_tr), by = 9L)] <- "c"
  fix$design$y_train <- factor(y_tr, levels = c("a", "b", "c"))
  fix$design$train_design$y <- fix$design$y_train

  y_te <- as.character(fix$design$y_test)
  y_te[seq(1L, length(y_te), by = 5L)] <- "z"
  fix$design$y_test <- factor(y_te, levels = c("a", "b", "c", "z"))
  fix$design$test_design$y <- fix$design$y_test

  old_opt <- options(
    rMVPA.searchlight_mode = "fast",
    rMVPA.naive_xdec_fast_kernel = NULL
  )
  on.exit(options(old_opt), add = TRUE)

  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = TRUE)
  out <- fit_roi(model, fix$roi_data, fix$context)

  expect_false(out$error)
  expect_true(all(is.finite(unname(out$result$probs)) | is.na(unname(out$result$probs))))
  expect_true(all(rowSums(out$result$probs, na.rm = FALSE) > 0 | is.na(rowSums(out$result$probs, na.rm = FALSE))))
  expect_equal(ncol(out$result$probs), length(levels(out$result$predicted)))
})

test_that("naive_xdec fast metrics match legacy performance computation (multiclass)", {
  fix <- build_naive_xdec_fast_fixture(seed = 9405L, nlevels = 3L)
  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = TRUE)

  core <- .naive_xdec_fit_core(
    Xtr = fix$roi_data$train_data,
    Xte = fix$roi_data$test_data,
    kernel = model$.fast_kernel
  )
  test_ind <- seq_len(nrow(fix$roi_data$test_data))

  fast_perf <- rMVPA:::.naive_xdec_fast_metrics(model, core, test_idx = test_ind)
  legacy_result <- classification_result(
    core$obs,
    core$pred,
    core$probs,
    testind = test_ind,
    test_design = fix$design$test_design,
    predictor = NULL
  )
  legacy_perf <- compute_performance(model, legacy_result)

  expect_true(is.numeric(fast_perf))
  expect_equal(fast_perf, legacy_perf, tolerance = 1e-12)
})

test_that("naive_xdec fast metrics match legacy performance computation (binary)", {
  fix <- build_naive_xdec_fast_fixture(seed = 9406L, nlevels = 2L)
  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = TRUE)

  core <- .naive_xdec_fit_core(
    Xtr = fix$roi_data$train_data,
    Xte = fix$roi_data$test_data,
    kernel = model$.fast_kernel
  )
  test_ind <- seq_len(nrow(fix$roi_data$test_data))

  fast_perf <- rMVPA:::.naive_xdec_fast_metrics(model, core, test_idx = test_ind)
  legacy_result <- classification_result(
    core$obs,
    core$pred,
    core$probs,
    testind = test_ind,
    test_design = fix$design$test_design,
    predictor = NULL
  )
  legacy_perf <- compute_performance(model, legacy_result)

  expect_true(is.numeric(fast_perf))
  expect_equal(fast_perf, legacy_perf, tolerance = 1e-12)
})

test_that("naive_xdec score-only metrics match legacy performance with splits", {
  split_by <- factor(rep(c("early", "late"), each = 18L))
  fix <- build_naive_xdec_fast_fixture(seed = 9409L, nlevels = 3L)
  fix$design$split_groups <- split(seq_len(length(split_by)), split_by)
  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = FALSE)

  core_scores <- .naive_xdec_fit_core(
    Xtr = fix$roi_data$train_data,
    Xte = fix$roi_data$test_data,
    kernel = model$.fast_kernel,
    return_probs = FALSE,
    return_scores = TRUE
  )
  test_ind <- seq_len(nrow(fix$roi_data$test_data))
  score_perf <- rMVPA:::.naive_xdec_fast_metrics_from_scores(model, core_scores, test_idx = test_ind)

  core_probs <- .naive_xdec_fit_core(
    Xtr = fix$roi_data$train_data,
    Xte = fix$roi_data$test_data,
    kernel = model$.fast_kernel
  )
  legacy_result <- classification_result(
    core_probs$obs,
    core_probs$pred,
    core_probs$probs,
    testind = test_ind,
    test_design = fix$design$test_design,
    predictor = NULL
  )
  legacy_perf <- compute_performance(model, legacy_result)

  expect_true(is.numeric(score_perf))
  expect_equal(score_perf, legacy_perf, tolerance = 1e-12)
})

test_that("naive_xdec score-only metrics match legacy performance for binary labels", {
  fix <- build_naive_xdec_fast_fixture(seed = 9410L, nlevels = 2L)
  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = FALSE)

  core_scores <- .naive_xdec_fit_core(
    Xtr = fix$roi_data$train_data,
    Xte = fix$roi_data$test_data,
    kernel = model$.fast_kernel,
    return_probs = FALSE,
    return_scores = TRUE
  )
  test_ind <- seq_len(nrow(fix$roi_data$test_data))
  score_perf <- rMVPA:::.naive_xdec_fast_metrics_from_scores(model, core_scores, test_idx = test_ind)

  core_probs <- .naive_xdec_fit_core(
    Xtr = fix$roi_data$train_data,
    Xte = fix$roi_data$test_data,
    kernel = model$.fast_kernel
  )
  legacy_result <- classification_result(
    core_probs$obs,
    core_probs$pred,
    core_probs$probs,
    testind = test_ind,
    test_design = fix$design$test_design,
    predictor = NULL
  )
  legacy_perf <- compute_performance(model, legacy_result)

  expect_true(is.numeric(score_perf))
  expect_equal(score_perf, legacy_perf, tolerance = 1e-12)
})

test_that("naive_xdec skips softmax when score-only fast metrics are available", {
  fix <- build_naive_xdec_fast_fixture(seed = 9411L, nlevels = 3L)
  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = FALSE)

  out <- testthat::with_mocked_bindings(
    fit_roi.naive_xdec_model(model, fix$roi_data, fix$context),
    .nx_softmax = function(...) stop(".nx_softmax should not be called"),
    .package = "rMVPA"
  )

  expect_false(out$error)
  expect_null(out$result)
  expect_true(is.numeric(out$metrics))
})

test_that("nx row correlation matches base cor and preserves p<2 degeneracy", {
  set.seed(9412L)
  Xa <- matrix(rnorm(8 * 5), nrow = 8)
  Xb <- matrix(rnorm(3 * 5), nrow = 3)

  expect_equal(.nx_row_cor(Xa, Xb), cor(t(Xa), t(Xb)), tolerance = 1e-12)

  deg <- .nx_row_cor(matrix(1:4, ncol = 1), matrix(5:7, ncol = 1))
  expect_true(all(is.na(deg)))
  expect_equal(dim(deg), c(4L, 3L))
})

test_that("naive_xdec skips classification_result when predictions are disabled", {
  fix <- build_naive_xdec_fast_fixture(seed = 9407L, nlevels = 3L)
  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = FALSE)

  out <- testthat::with_mocked_bindings(
    fit_roi.naive_xdec_model(model, fix$roi_data, fix$context),
    classification_result = function(...) stop("classification_result should not be called"),
    .package = "rMVPA"
  )

  expect_false(out$error)
  expect_null(out$result)
  expect_true(is.numeric(out$metrics))
})

test_that("naive_xdec falls back to legacy performance path for custom metrics", {
  fix <- build_naive_xdec_fast_fixture(seed = 9408L, nlevels = 3L)
  model <- naive_xdec_model(fix$dataset, fix$design, return_predictions = FALSE)
  model$performance <- function(result) {
    c(custom_acc = mean(result$observed == result$predicted))
  }

  called <- FALSE
  original_classification_result <- getFromNamespace("classification_result", "rMVPA")

  out <- testthat::with_mocked_bindings(
    fit_roi.naive_xdec_model(model, fix$roi_data, fix$context),
    classification_result = function(...) {
      called <<- TRUE
      original_classification_result(...)
    },
    .package = "rMVPA"
  )

  expect_true(called)
  expect_false(out$error)
  expect_true("custom_acc" %in% names(out$metrics))
})
