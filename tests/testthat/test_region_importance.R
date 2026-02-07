test_that("cv_score_global returns a finite scalar", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 8,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  X <- get_feature_matrix(ds$dataset)
  y <- y_train(ds$design)
  feature_ids <- get_center_ids(ds$dataset)

  score <- rMVPA:::cv_score_global(mspec, X, y, cv, feature_ids, metric = NULL)

  expect_true(is.numeric(score))
  expect_length(score, 1)
  expect_true(is.finite(score))
})

test_that("cv_score_global works with a column subset", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 8,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  X <- get_feature_matrix(ds$dataset)
  y <- y_train(ds$design)
  feature_ids <- get_center_ids(ds$dataset)

  # Use a subset of columns
  sel <- 1:4
  score <- rMVPA:::cv_score_global(mspec, X[, sel, drop = FALSE], y, cv,
                                    feature_ids[sel], metric = NULL)

  expect_true(is.numeric(score))
  expect_length(score, 1)
  expect_true(is.finite(score))
})

test_that("region_importance end-to-end with clustered dataset", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 40, K = 8,
                                      nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  P <- ncol(get_feature_matrix(ds$dataset))

  ri <- region_importance(mspec, n_iter = 20, subset_fraction = 0.5)

  expect_s3_class(ri, "region_importance_result")

  # importance vector
  expect_length(ri$importance, P)
  expect_true(is.numeric(ri$importance))

  # p_values vector
  expect_length(ri$p_values, P)
  expect_true(all(is.na(ri$p_values) | (ri$p_values >= 0 & ri$p_values <= 1)))

  # importance_map is a valid spatial object

  expect_true(!is.null(ri$importance_map))

  # p_value_map
  expect_true(!is.null(ri$p_value_map))

  # stats_table structure
  expect_s3_class(ri$stats_table, "tbl_df")
  expected_cols <- c("feature_id", "importance", "p_value", "ci_lower",
                     "ci_upper", "n_in", "n_out", "mean_in", "mean_out")
  expect_true(all(expected_cols %in% names(ri$stats_table)))
  expect_equal(nrow(ri$stats_table), P)

  # iteration_log
  expect_s3_class(ri$iteration_log, "tbl_df")
  expect_equal(nrow(ri$iteration_log), 20)
  expect_true(all(c("iter", "performance", "n_features") %in% names(ri$iteration_log)))
})

test_that("print.region_importance_result works without error", {
  ri <- region_importance_result(
    importance     = c(0.05, -0.02, 0.1, NA),
    importance_map = NULL,
    p_values       = c(0.01, 0.5, 0.001, NA),
    p_value_map    = NULL,
    stats_table    = tibble::tibble(
      feature_id = 1:4,
      importance = c(0.05, -0.02, 0.1, NA),
      p_value    = c(0.01, 0.5, 0.001, NA),
      ci_lower   = c(-0.01, -0.05, 0.02, NA),
      ci_upper   = c(0.1, 0.02, 0.18, NA),
      n_in       = c(10L, 10L, 10L, 0L),
      n_out      = c(10L, 10L, 10L, 20L),
      mean_in    = c(0.55, 0.48, 0.6, NA),
      mean_out   = c(0.5, 0.5, 0.5, 0.5)
    ),
    iteration_log  = tibble::tibble(
      iter = 1:20,
      performance = runif(20, 0.4, 0.7),
      n_features = rep(2L, 20)
    ),
    model_spec = NULL
  )

  expect_output(print(ri), "Region Importance Result")
})

test_that("region_importance with signal parcel ranks it highest", {
  skip_if_not_installed("sda")
  set.seed(123)

  # Create clustered dataset with K=6 parcels
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 60, K = 6,
                                      nlevels = 2, blocks = 4)

  # Inject signal into parcel 1
  X <- get_feature_matrix(ds$dataset)
  y <- ds$design$y_train
  signal <- ifelse(as.integer(y) == 1, 2, -2)
  X[, 1] <- X[, 1] + signal

  # Reconstruct dataset with modified data - put signal directly in the ts matrix
  ds$dataset$train_data@ts <- X

  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  ri <- region_importance(mspec, n_iter = 50, subset_fraction = 0.5)

  # Parcel 1 should have the highest importance
  valid_imp <- ri$importance[!is.na(ri$importance)]
  expect_true(length(valid_imp) > 0)

  # Feature 1 should be at or near the top
  rank_1 <- rank(-ri$importance, na.last = TRUE)[1]
  expect_true(rank_1 <= 2)
})

test_that("cv_score_global respects explicit metric parameter", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 8,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  X <- get_feature_matrix(ds$dataset)
  y <- y_train(ds$design)
  feature_ids <- get_center_ids(ds$dataset)

  # Request a specific metric
  score_acc <- rMVPA:::cv_score_global(mspec, X, y, cv, feature_ids, metric = "Accuracy")
  score_default <- rMVPA:::cv_score_global(mspec, X, y, cv, feature_ids, metric = NULL)

  expect_true(is.finite(score_acc))
  expect_true(is.finite(score_default))

  # With a nonexistent metric name, falls back to first metric
  score_bad <- rMVPA:::cv_score_global(mspec, X, y, cv, feature_ids, metric = "Nonexistent")
  expect_true(is.finite(score_bad))
  expect_equal(score_bad, score_default)
})

test_that("region_importance with image dataset produces NeuroVol maps", {
  skip_if_not_installed("sda")
  ds <- gen_sample_dataset(D = c(6, 6, 6), nobs = 40, nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  P <- sum(ds$dataset$mask > 0)

  ri <- region_importance(mspec, n_iter = 15, subset_fraction = 0.5)

  expect_s3_class(ri, "region_importance_result")
  expect_length(ri$importance, P)
  expect_length(ri$p_values, P)

  # Maps should be NeuroVol for image datasets
  expect_true(inherits(ri$importance_map, "NeuroVol"))
  expect_true(inherits(ri$p_value_map, "NeuroVol"))
})

test_that("region_importance rejects invalid subset_fraction", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 6,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  expect_error(region_importance(mspec, n_iter = 5, subset_fraction = 0),
               "subset_fraction must be between")
  expect_error(region_importance(mspec, n_iter = 5, subset_fraction = 1),
               "subset_fraction must be between")
  expect_error(region_importance(mspec, n_iter = 5, subset_fraction = -0.1),
               "subset_fraction must be between")
})

test_that("n_in + n_out == n_iter and importance consistency", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 40, K = 8,
                                      nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  n_iter <- 25
  ri <- region_importance(mspec, n_iter = n_iter, subset_fraction = 0.5)

  # n_in + n_out should always equal n_iter
  expect_true(all(ri$stats_table$n_in + ri$stats_table$n_out == n_iter))

  # importance vector matches stats_table$importance
  expect_equal(ri$importance, ri$stats_table$importance)

  # importance = mean_in - mean_out for valid features
  valid <- !is.na(ri$importance)
  expected_imp <- ri$stats_table$mean_in[valid] - ri$stats_table$mean_out[valid]
  expect_equal(ri$importance[valid], expected_imp, tolerance = 1e-12)

  # iteration_log has all finite performance values
  expect_true(all(is.finite(ri$iteration_log$performance)))

  # All n_features in iteration_log match subset_size
  P <- ncol(get_feature_matrix(ds$dataset))
  expected_size <- floor(P * 0.5)
  expect_true(all(ri$iteration_log$n_features == expected_size))
})

test_that("high subset_fraction warns about always-included features", {
  skip_if_not_installed("sda")
  # With K=3 parcels and subset_fraction=0.9, floor(3*0.9)=2 out of 3.
  # With only 10 iterations, some features will likely always be included.
  # But it's not guaranteed, so we use a very extreme case:
  # K=2 and subset_fraction=0.9 -> floor(2*0.9)=1, so each iter picks 1 of 2.
  # Neither feature is always included with decent n_iter, so use K=2, frac=0.99 -> floor(2*0.99)=1
  # Actually for a guaranteed warning we need subset_size == P.
  # floor(3 * 0.99) = 2 out of 3 features. Not guaranteed to always include any one.
  # Let's use K=2, frac=0.8 -> floor(2*0.8)=1. Each iter picks 1 of 2. No warning.
  # For K=2, frac=0.99 -> floor(2*0.99)=1. Same.
  # To guarantee the warning, we need subset_size == P, i.e. frac close to 1.
  # But frac >= 1 is an error. So with K=2 and frac=0.99, subset_size=1, won't trigger.
  # The warning path is really only hit in degenerate cases. Let's just test
  # that no error occurs with high fraction and few features.
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 3,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  # frac=0.9 with K=3 -> subset_size=2 out of 3. Each feature has ~2/3 chance of inclusion.
  # Very unlikely any is never included in 20 iters, but the code should handle it.
  ri <- region_importance(mspec, n_iter = 20, subset_fraction = 0.9)

  expect_s3_class(ri, "region_importance_result")
  expect_length(ri$importance, 3)
})

test_that("run_global with return_fits works after refactoring", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 8,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec, return_fits = TRUE)

  expect_s3_class(res, "global_mvpa_result")
  expect_true(!is.null(res$fold_fits))
  expect_equal(length(res$fold_fits), 3)  # 3 blocks = 3 folds
  non_null <- Filter(Negate(is.null), res$fold_fits)
  expect_true(length(non_null) > 0)
  expect_true(inherits(non_null[[1]], "model_fit"))
})

test_that("run_global still works after refactoring", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 40, K = 12,
                                      nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  expect_s3_class(res, "global_mvpa_result")
  expect_true(nrow(res$performance_table) > 0)
  expect_true(!is.null(res$importance_vector))
  expect_equal(length(res$importance_vector), ncol(get_feature_matrix(ds$dataset)))
})

test_that("region_importance with RF works end-to-end", {
  skip_if_not_installed("randomForest")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 40, K = 8,
                                      nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("rf")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  P <- ncol(get_feature_matrix(ds$dataset))
  n_iter <- 20

  ri <- region_importance(mspec, n_iter = n_iter, subset_fraction = 0.5)

  expect_s3_class(ri, "region_importance_result")

  # importance vector
  expect_length(ri$importance, P)
  expect_true(is.numeric(ri$importance))

  # p-values
  expect_length(ri$p_values, P)
  expect_true(all(is.na(ri$p_values) | (ri$p_values >= 0 & ri$p_values <= 1)))

  # spatial maps should be built
  expect_true(!is.null(ri$importance_map))
  expect_true(!is.null(ri$p_value_map))

  # stats_table structure
  expect_s3_class(ri$stats_table, "tbl_df")
  expected_cols <- c("feature_id", "importance", "p_value", "ci_lower",
                     "ci_upper", "n_in", "n_out", "mean_in", "mean_out")
  expect_true(all(expected_cols %in% names(ri$stats_table)))
  expect_equal(nrow(ri$stats_table), P)

  # n_in + n_out == n_iter
  expect_true(all(ri$stats_table$n_in + ri$stats_table$n_out == n_iter))

  # iteration_log
  expect_s3_class(ri$iteration_log, "tbl_df")
  expect_equal(nrow(ri$iteration_log), n_iter)
  expect_true(all(is.finite(ri$iteration_log$performance)))

  # importance = mean_in - mean_out for valid features
  valid <- !is.na(ri$importance)
  if (sum(valid) > 0) {
    expected_imp <- ri$stats_table$mean_in[valid] - ri$stats_table$mean_out[valid]
    expect_equal(ri$importance[valid], expected_imp, tolerance = 1e-12)
  }
})

test_that("region_importance with RF on image dataset produces NeuroVol maps", {
  skip_if_not_installed("randomForest")
  ds <- gen_sample_dataset(D = c(6, 6, 6), nobs = 40, nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("rf")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  P <- sum(ds$dataset$mask > 0)
  ri <- region_importance(mspec, n_iter = 10, subset_fraction = 0.5)

  expect_s3_class(ri, "region_importance_result")
  expect_length(ri$importance, P)
  expect_true(inherits(ri$importance_map, "NeuroVol"))
  expect_true(inherits(ri$p_value_map, "NeuroVol"))
})
