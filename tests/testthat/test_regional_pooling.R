context("regional pooling and stacking")

make_classification_predtab <- function(n = 50, nroi = 4, seed = 1, rownum_offset = 0L) {
  set.seed(seed)
  rownum <- as.integer(rownum_offset + seq_len(n))
  observed <- factor(sample(c("A", "B"), n, replace = TRUE), levels = c("A", "B"))
  signal <- rnorm(n)

  tabs <- lapply(seq_len(nroi), function(r) {
    eta <- signal + rnorm(n, sd = 0.8) + (r - 1) * 0.1
    prob_b <- plogis(eta)
    tibble::tibble(
      .rownum = rownum,
      roinum = as.integer(r),
      observed = observed,
      prob_A = 1 - prob_b,
      prob_B = prob_b,
      predicted = ifelse(prob_b > 0.5, "B", "A"),
      correct = predicted == as.character(observed)
    )
  })

  dplyr::bind_rows(tabs)
}

make_regression_predtab <- function(n = 60, nroi = 5, seed = 2, rownum_offset = 0L) {
  set.seed(seed)
  rownum <- as.integer(rownum_offset + seq_len(n))
  observed <- rnorm(n)
  signal <- observed + rnorm(n, sd = 0.5)

  tabs <- lapply(seq_len(nroi), function(r) {
    pred <- signal + rnorm(n, sd = 0.6) + (r - 1) * 0.05
    tibble::tibble(
      .rownum = rownum,
      roinum = as.integer(r),
      observed = observed,
      predicted = pred
    )
  })

  dplyr::bind_rows(tabs)
}

test_that("weighted mean pooling uses provided weights", {
  observed <- factor(c("A", "B", "A"), levels = c("A", "B"))
  p1 <- tibble::tibble(
    .rownum = 1:3,
    roinum = 1L,
    observed = observed,
    prob_A = c(0.9, 0.1, 0.8),
    prob_B = c(0.1, 0.9, 0.2),
    predicted = c("A", "B", "A"),
    correct = TRUE
  )
  p2 <- tibble::tibble(
    .rownum = 1:3,
    roinum = 2L,
    observed = observed,
    prob_A = c(0.3, 0.7, 0.2),
    prob_B = c(0.7, 0.3, 0.8),
    predicted = c("B", "A", "B"),
    correct = FALSE
  )

  pooled <- combine_prediction_tables(list(p1, p2), wts = c(3, 1), collapse_regions = TRUE)
  expected_prob_a <- (3 * p1$prob_A + p2$prob_A) / 4

  expect_equal(pooled$prob_A, expected_prob_a, tolerance = 1e-12)
  expect_equal(pooled$prob_B, 1 - expected_prob_a, tolerance = 1e-12)
})

test_that("stacking is deterministic with stack_seed", {
  stack_fun <- getFromNamespace(".crossfit_stack_predictions", "rMVPA")
  ptab <- make_classification_predtab(n = 40, nroi = 3, seed = 11)
  model_spec <- list(design = list(block_var = NULL))

  s1 <- stack_fun(ptab, model_spec, stack_folds = 5, stack_seed = 101)
  s2 <- stack_fun(ptab, model_spec, stack_folds = 5, stack_seed = 101)
  s3 <- stack_fun(ptab, model_spec, stack_folds = 5, stack_seed = 202)

  expect_equal(s1, s2, tolerance = 1e-12)
  expect_true(any(abs(s1$prob_A - s3$prob_A) > 1e-8))
})

test_that("stacking supports integer, fold-id vector, and fold-list specs", {
  stack_fun <- getFromNamespace(".crossfit_stack_predictions", "rMVPA")
  ptab <- make_classification_predtab(n = 36, nroi = 3, seed = 33, rownum_offset = 100L)
  model_spec <- list(design = list(block_var = NULL))

  rownum_unique <- sort(unique(ptab$.rownum))
  fold_ids <- rep(1:6, length.out = length(rownum_unique))
  fold_list_idx <- split(seq_along(rownum_unique), fold_ids)
  fold_list_rownum <- lapply(fold_list_idx, function(i) rownum_unique[i])

  s_vec <- stack_fun(ptab, model_spec, stack_folds = fold_ids)
  s_idx <- stack_fun(ptab, model_spec, stack_folds = fold_list_idx)
  s_row <- stack_fun(ptab, model_spec, stack_folds = fold_list_rownum)

  expect_equal(s_vec$prob_A, s_idx$prob_A, tolerance = 1e-12)
  expect_equal(s_vec$prob_A, s_row$prob_A, tolerance = 1e-12)
  expect_equal(s_vec$.rownum, rownum_unique)
})

test_that("stacking for regression returns expected columns and metrics", {
  stack_fun <- getFromNamespace(".crossfit_stack_predictions", "rMVPA")
  perf_fun <- getFromNamespace(".compute_pooled_performance", "rMVPA")
  ptab <- make_regression_predtab(n = 45, nroi = 4, seed = 44)
  model_spec <- list(design = list(block_var = NULL))

  pooled <- stack_fun(ptab, model_spec, stack_folds = 5, stack_seed = 99)
  perf <- perf_fun(pooled)

  expect_true(all(c(".rownum", "roinum", "observed", "predicted", "residual", "abs_error", "sq_error") %in% names(pooled)))
  expect_true(all(c("R2", "RMSE", "spearcor") %in% names(perf)))
})

test_that("pool_predictions mean validates pooled_weights length", {
  pool_fun <- getFromNamespace(".pool_regional_predictions", "rMVPA")
  ptab <- make_classification_predtab(n = 20, nroi = 3, seed = 55)
  model_spec <- list(design = list(block_var = NULL))

  expect_error(
    pool_fun(ptab, model_spec, method = "mean", pooled_weights = c(1, 2)),
    "`pooled_weights` must have one weight per ROI"
  )
})

test_that("cross-fitted stacking does not trivially overfit random labels", {
  stack_fun <- getFromNamespace(".crossfit_stack_predictions", "rMVPA")
  model_spec <- list(design = list(block_var = NULL))
  ptab <- make_classification_predtab(n = 80, nroi = 6, seed = 88)

  # Replace observed labels with balanced random labels unrelated to ROI predictions.
  set.seed(123)
  rownum_unique <- sort(unique(ptab$.rownum))
  y_rand <- sample(rep(c("A", "B"), each = length(rownum_unique) / 2))
  map_y <- stats::setNames(y_rand, rownum_unique)
  ptab$observed <- factor(map_y[as.character(ptab$.rownum)], levels = c("A", "B"))
  ptab$correct <- ptab$predicted == as.character(ptab$observed)

  s <- stack_fun(ptab, model_spec, stack_folds = 5, stack_seed = 42)
  acc <- mean(s$predicted == as.character(s$observed))

  # On random labels, OOF stacking should not show near-perfect performance.
  expect_lt(acc, 0.8)
})
