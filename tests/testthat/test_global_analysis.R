test_that("haufe_importance computes correct activation patterns", {
  set.seed(1)
  P <- 5
  D <- 2
  W <- matrix(rnorm(P * D), P, D)
  Sigma_x <- crossprod(matrix(rnorm(20 * P), 20, P)) / 20  # PSD matrix

  result <- haufe_importance(W, Sigma_x)

  # A = Sigma_x %*% W %*% solve(t(W) %*% Sigma_x %*% W)
  WtSW <- t(W) %*% Sigma_x %*% W
  expected_A <- Sigma_x %*% W %*% solve(WtSW)

  expect_equal(result$A, expected_A, tolerance = 1e-10)
  expect_equal(length(result$importance), P)
  expect_true(all(result$importance >= 0))
  # L2 norm check
  expected_imp <- sqrt(rowSums(expected_A^2))
  expect_equal(result$importance, expected_imp, tolerance = 1e-10)
})

test_that("haufe_importance uses ginv for singular matrices", {
  P <- 3
  D <- 2
  # Make W such that W'*Sigma_x*W is singular (W columns are linearly dependent)
  W <- matrix(c(1, 0, 0, 2, 0, 0), P, D)  # col2 = 2*col1
  Sigma_x <- diag(P)

  # Should not error (falls back to MASS::ginv)
  result <- haufe_importance(W, Sigma_x)
  expect_equal(nrow(result$A), P)
  expect_equal(ncol(result$A), D)
  expect_equal(length(result$importance), P)
})

test_that("haufe_importance accepts custom summary_fun", {
  set.seed(2)
  P <- 4
  D <- 2
  W <- matrix(rnorm(P * D), P, D)
  Sigma_x <- diag(P)

  # Use max absolute value instead of L2
  result <- haufe_importance(W, Sigma_x,
                              summary_fun = function(A) apply(abs(A), 1, max))
  expect_equal(length(result$importance), P)
})

test_that("extract_weights.sda returns correct dimensions", {
  skip_if_not_installed("sda")
  set.seed(10)
  X <- matrix(rnorm(50 * 10), 50, 10)
  y <- factor(rep(letters[1:3], length.out = 50))

  fit <- sda::sda(X, y, verbose = FALSE)
  W <- extract_weights(fit)

  expect_true(is.matrix(W))
  expect_equal(nrow(W), 10)      # P features
  expect_equal(ncol(W), 3)       # K classes for sda beta
})

test_that("extract_weights.glmnet works for binary classification", {
  skip_if_not_installed("glmnet")
  set.seed(20)
  X <- matrix(rnorm(100 * 8), 100, 8)
  y <- factor(rep(c("a", "b"), 50))

  fit <- glmnet::glmnet(X, y, family = "binomial", nlambda = 10)
  fit$opt_lambda <- fit$lambda[5]

  W <- extract_weights(fit)
  expect_true(is.matrix(W))
  expect_equal(nrow(W), 8)
  expect_equal(ncol(W), 1)
})

test_that("extract_weights.glmnet works for multinomial classification", {
  skip_if_not_installed("glmnet")
  set.seed(30)
  X <- matrix(rnorm(90 * 6), 90, 6)
  y <- factor(rep(c("a", "b", "c"), 30))

  fit <- glmnet::glmnet(X, y, family = "multinomial", nlambda = 10)
  fit$opt_lambda <- fit$lambda[5]

  W <- extract_weights(fit)
  expect_true(is.matrix(W))
  expect_equal(nrow(W), 6)
  expect_equal(ncol(W), 3)   # one column per class for multinomial
})

test_that("extract_weights.default errors informatively", {
  obj <- list(foo = 1)
  class(obj) <- "weird_model"
  expect_error(extract_weights(obj), "no method for class")
})

test_that("get_feature_matrix.matrix is identity", {
  M <- matrix(1:12, 3, 4)
  expect_identical(get_feature_matrix(M), M)
})

test_that("get_feature_matrix.mvpa_image_dataset returns correct dims", {
  ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 10, blocks = 2)
  X <- get_feature_matrix(ds$dataset)
  expect_true(is.matrix(X))
  expect_equal(nrow(X), 10)
  expect_equal(ncol(X), sum(ds$dataset$mask > 0))
})

test_that("get_feature_matrix.mvpa_clustered_dataset returns correct dims", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 8, blocks = 3)
  X <- get_feature_matrix(ds$dataset)
  expect_true(is.matrix(X))
  expect_equal(nrow(X), 20)
  expect_equal(ncol(X), neuroim2::num_clusters(ds$dataset$cvol))
})

test_that("run_global with clustered dataset + sda_notune", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 40, K = 12,
                                      nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  expect_s3_class(res, "global_mvpa_result")
  expect_true(nrow(res$performance_table) > 0)
  expect_true(!is.null(res$result))
  expect_true(!is.null(res$importance_vector))
  expect_equal(length(res$importance_vector), ncol(get_feature_matrix(ds$dataset)))
  expect_true(!is.null(res$importance_map))
  expect_true(!is.null(res$activation_patterns))
  expect_true(!is.null(res$raw_weights))
})

test_that("run_global with image dataset + sda_notune", {
  skip_if_not_installed("sda")
  ds <- gen_sample_dataset(D = c(6, 6, 6), nobs = 40, nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  expect_s3_class(res, "global_mvpa_result")
  expect_true(nrow(res$performance_table) > 0)
  P <- sum(ds$dataset$mask > 0)
  expect_equal(length(res$importance_vector), P)
  expect_true(inherits(res$importance_map, "NeuroVol"))
})

test_that("run_global with plain matrix X", {
  skip_if_not_installed("sda")
  ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 30, nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  # Use a smaller subset of features as plain matrix
  X <- get_feature_matrix(ds$dataset)[, 1:10]

  res <- run_global(mspec, X = X)

  expect_s3_class(res, "global_mvpa_result")
  expect_equal(length(res$importance_vector), 10)
})

test_that("run_global return_fits stores per-fold fits", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 8,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec, return_fits = TRUE)

  expect_true(!is.null(res$fold_fits))
  expect_equal(length(res$fold_fits), 3)  # 3 blocks = 3 folds
  # Each fold fit should be a model_fit object
  non_null <- Filter(Negate(is.null), res$fold_fits)
  expect_true(length(non_null) > 0)
  expect_true(inherits(non_null[[1]], "model_fit"))
})

test_that("feature mask zero-padding works correctly", {
  skip_if_not_installed("sda")
  # Use dataset with some NA columns to force feature_mask filtering
  ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 30, nlevels = 2,
                           blocks = 3, na_cols = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  # Importance vector should have length == total features (including NAs)
  P <- sum(ds$dataset$mask > 0)
  expect_equal(length(res$importance_vector), P)
})

test_that("print.global_mvpa_result works without error", {
  res <- global_mvpa_result(
    performance_table   = tibble::tibble(Accuracy = 0.75),
    result              = NULL,
    importance_map      = NULL,
    importance_vector   = c(0.1, 0.2, 0.3),
    activation_patterns = matrix(rnorm(6), 3, 2),
    raw_weights         = matrix(rnorm(6), 3, 2),
    fold_fits           = NULL,
    model_spec          = NULL
  )

  expect_output(print(res), "Global MVPA Result")
})

# ===== model_importance unit tests =====

test_that("model_importance.sda matches haufe_importance directly", {
  skip_if_not_installed("sda")
  set.seed(42)
  P <- 8
  X <- matrix(rnorm(50 * P), 50, P)
  y <- factor(rep(letters[1:2], 25))
  fit <- sda::sda(X, y, verbose = FALSE)

  # model_importance should give same result as manual haufe_importance
  imp <- model_importance(fit, X)
  W <- extract_weights(fit)
  Sigma_x <- cov(X)
  expected <- haufe_importance(W, Sigma_x)

  expect_true(is.numeric(imp))
  expect_equal(length(imp), P)
  expect_equal(imp, expected$importance, tolerance = 1e-10)
  expect_true(all(imp >= 0))  # L2 norm is non-negative
})

test_that("model_importance.sda with custom summary_fun", {
  skip_if_not_installed("sda")
  set.seed(42)
  X <- matrix(rnorm(50 * 6), 50, 6)
  y <- factor(rep(letters[1:3], length.out = 50))
  fit <- sda::sda(X, y, verbose = FALSE)

  custom_fun <- function(A) apply(abs(A), 1, max)
  imp_custom <- model_importance(fit, X, summary_fun = custom_fun)
  imp_default <- model_importance(fit, X)

  expect_equal(length(imp_custom), 6)
  # Custom and default should generally differ
  expect_false(isTRUE(all.equal(imp_custom, imp_default)))
})

test_that("model_importance.glmnet matches haufe_importance directly", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  P <- 10
  X <- matrix(rnorm(100 * P), 100, P)
  y <- factor(rep(c("a", "b"), 50))

  fit <- glmnet::glmnet(X, y, family = "binomial", nlambda = 10)
  fit$opt_lambda <- fit$lambda[5]

  imp <- model_importance(fit, X)
  W <- extract_weights(fit)
  Sigma_x <- cov(X)
  expected <- haufe_importance(W, Sigma_x)

  expect_true(is.numeric(imp))
  expect_equal(length(imp), P)
  expect_equal(imp, expected$importance, tolerance = 1e-10)
})

test_that("model_importance.glmnet works for multinomial", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  P <- 8
  X <- matrix(rnorm(90 * P), 90, P)
  y <- factor(rep(c("a", "b", "c"), 30))

  fit <- glmnet::glmnet(X, y, family = "multinomial", nlambda = 10)
  fit$opt_lambda <- fit$lambda[5]

  imp <- model_importance(fit, X)
  expect_true(is.numeric(imp))
  expect_equal(length(imp), P)
  expect_true(all(imp >= 0))
})

test_that("model_importance.randomForest with importance=TRUE uses MeanDecreaseAccuracy", {
  skip_if_not_installed("randomForest")
  set.seed(42)
  P <- 5
  X <- matrix(rnorm(60 * P), 60, P)
  y <- factor(rep(letters[1:3], 20))
  fit <- randomForest::randomForest(X, y, importance = TRUE)

  imp <- model_importance(fit, X)
  raw_imp <- randomForest::importance(fit)
  expected <- as.numeric(raw_imp[, "MeanDecreaseAccuracy"])

  expect_true(is.numeric(imp))
  expect_equal(length(imp), P)
  expect_equal(imp, expected)
})

test_that("model_importance.randomForest without importance=TRUE uses MeanDecreaseGini", {
  skip_if_not_installed("randomForest")
  set.seed(42)
  P <- 5
  X <- matrix(rnorm(60 * P), 60, P)
  y <- factor(rep(letters[1:3], 20))
  # Default: importance=FALSE, so only Gini is available
  fit <- randomForest::randomForest(X, y)

  imp <- model_importance(fit, X)
  raw_imp <- randomForest::importance(fit)
  expected <- as.numeric(raw_imp[, "MeanDecreaseGini"])

  expect_true(is.numeric(imp))
  expect_equal(length(imp), P)
  expect_equal(imp, expected)
})

test_that("model_importance.default returns NULL", {
  obj <- list(foo = 1)
  class(obj) <- "unknown_model"
  expect_null(model_importance(obj, matrix(0, 2, 2)))
})

# ===== run_global + model_importance integration tests =====

test_that("run_global with RF produces importance map (clustered dataset)", {
  skip_if_not_installed("randomForest")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 40, K = 12,
                                      nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("rf")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  expect_s3_class(res, "global_mvpa_result")
  expect_true(nrow(res$performance_table) > 0)
  # RF importance_vector should be non-NULL, correct length, all positive (Gini)
  expect_true(!is.null(res$importance_vector))
  expect_equal(length(res$importance_vector), ncol(get_feature_matrix(ds$dataset)))
  expect_true(all(res$importance_vector >= 0))
  expect_true(!is.null(res$importance_map))
  # RF has no linear weights or activation patterns
  expect_true(is.null(res$raw_weights))
  expect_true(is.null(res$activation_patterns))
})

test_that("run_global with RF produces importance map (image dataset)", {
  skip_if_not_installed("randomForest")
  ds <- gen_sample_dataset(D = c(6, 6, 6), nobs = 40, nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("rf")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  expect_s3_class(res, "global_mvpa_result")
  P <- sum(ds$dataset$mask > 0)
  expect_true(!is.null(res$importance_vector))
  expect_equal(length(res$importance_vector), P)
  expect_true(inherits(res$importance_map, "NeuroVol"))
  expect_true(is.null(res$raw_weights))
  expect_true(is.null(res$activation_patterns))
})

test_that("run_global with RF + return_fits stores per-fold fits", {
  skip_if_not_installed("randomForest")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 8,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("rf")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec, return_fits = TRUE)

  expect_s3_class(res, "global_mvpa_result")
  expect_true(!is.null(res$fold_fits))
  expect_equal(length(res$fold_fits), 3)
  non_null <- Filter(Negate(is.null), res$fold_fits)
  expect_true(length(non_null) > 0)
  expect_true(inherits(non_null[[1]], "model_fit"))
  # Importance should still be computed alongside return_fits
  expect_true(!is.null(res$importance_vector))
})

test_that("feature mask zero-padding works for model_importance path", {
  skip_if_not_installed("sda")
  # Dataset with NA columns triggers feature masking in train_model
  ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 30, nlevels = 2,
                           blocks = 3, na_cols = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  P <- sum(ds$dataset$mask > 0)
  expect_equal(length(res$importance_vector), P)
  # Masked-out features should have zero importance (zero-padded)
  expect_true(any(res$importance_vector == 0))
  # Non-masked features should have non-zero importance
  expect_true(any(res$importance_vector != 0))
})

test_that("SDA run_global still produces activation_patterns and raw_weights", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 40, K = 10,
                                      nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec)

  P <- ncol(get_feature_matrix(ds$dataset))
  # importance_vector from per-fold Haufe
  expect_true(!is.null(res$importance_vector))
  expect_equal(length(res$importance_vector), P)
  # backward-compat: raw_weights and activation_patterns still populated
  expect_true(!is.null(res$raw_weights))
  expect_true(is.matrix(res$raw_weights))
  expect_equal(nrow(res$raw_weights), P)
  expect_true(!is.null(res$activation_patterns))
  expect_true(is.matrix(res$activation_patterns))
  expect_equal(nrow(res$activation_patterns), P)
})

test_that("cv_run_global compute_importance=FALSE returns empty fold_importance", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 6,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  X <- get_feature_matrix(ds$dataset)
  y <- y_train(ds$design)
  fids <- get_center_ids(ds$dataset)

  cv_out <- rMVPA:::cv_run_global(mspec, X, y, cv, fids,
                                   compute_importance = FALSE)

  # All fold_importance entries should be NULL when compute_importance=FALSE
  expect_true(all(sapply(cv_out$fold_importance, is.null)))
})

test_that("cv_run_global compute_importance=TRUE returns per-fold vectors", {
  skip_if_not_installed("sda")
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 6,
                                      nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("sda_notune")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  X <- get_feature_matrix(ds$dataset)
  P <- ncol(X)
  y <- y_train(ds$design)
  fids <- get_center_ids(ds$dataset)

  cv_out <- rMVPA:::cv_run_global(mspec, X, y, cv, fids,
                                   compute_importance = TRUE)

  valid_imp <- Filter(Negate(is.null), cv_out$fold_importance)
  expect_true(length(valid_imp) > 0)
  # Each per-fold importance vector should be P-length
  for (imp_k in valid_imp) {
    expect_equal(length(imp_k), P)
    expect_true(is.numeric(imp_k))
  }
})
