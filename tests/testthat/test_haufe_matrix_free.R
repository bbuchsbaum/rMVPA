library(testthat)

test_that("matrix-free Haufe equals the dense covariance reference", {
  set.seed(11)
  n <- 60; P <- 25; D <- 3
  X <- matrix(rnorm(n * P), n, P)
  X[, 1:5] <- X[, 1:5] + 3   # non-zero means so centring matters
  W <- matrix(rnorm(P * D), P, D)

  dense <- haufe_importance(W, cov(X))
  free  <- haufe_importance(W, X = X)

  expect_equal(unname(free$A), unname(dense$A), tolerance = 1e-10)
  expect_equal(free$importance, dense$importance, tolerance = 1e-10)
})

test_that("matrix-free Haufe blockwise accumulation matches the single pass", {
  set.seed(12)
  X <- matrix(rnorm(40 * 33), 40, 33)
  W <- matrix(rnorm(33 * 2), 33, 2)
  one   <- haufe_importance(W, X = X)
  block <- haufe_importance(W, X = X, block_size = 7)
  expect_equal(block$A, one$A, tolerance = 1e-12)
  expect_equal(block$importance, one$importance, tolerance = 1e-12)
})

test_that("matrix-free Haufe handles rank-deficient W like the dense path", {
  set.seed(13)
  X <- matrix(rnorm(50 * 12), 50, 12)
  w <- rnorm(12)
  W <- cbind(w, 2 * w)   # rank 1
  dense <- haufe_importance(W, cov(X))
  free  <- haufe_importance(W, X = X)
  expect_equal(unname(free$A), unname(dense$A), tolerance = 1e-8)
})

test_that("matrix-free Haufe honours center = FALSE against an uncentred cross-product", {
  set.seed(14)
  X <- matrix(rnorm(30 * 8), 30, 8) + 1
  W <- matrix(rnorm(8 * 2), 8, 2)
  S_unc <- crossprod(X) / (nrow(X) - 1)
  dense <- haufe_importance(W, S_unc)
  free  <- haufe_importance(W, X = X, center = FALSE)
  expect_equal(unname(free$A), unname(dense$A), tolerance = 1e-10)
})

test_that("haufe_importance validates its inputs", {
  W <- matrix(rnorm(6), 3, 2)
  expect_error(haufe_importance(W), "Sigma_x")
  expect_error(haufe_importance(W, X = matrix(rnorm(20), 5, 4)), "nrow\\(W\\)")
})

test_that("model_importance methods no longer allocate a P x P covariance", {
  skip_if_not_installed("sda")
  set.seed(15)
  n <- 30; P <- 400
  X <- matrix(rnorm(n * P), n, P)
  y <- factor(rep(c("a", "b"), length.out = n))
  fit <- sda::sda(X, y, verbose = FALSE)
  # Structural guard: none of the Haufe-based importance methods, nor the
  # global summary path, may form a dense feature covariance.
  bodies <- vapply(
    list(rMVPA:::model_importance.sda, rMVPA:::model_importance.glmnet,
         rMVPA:::model_importance.spacenet_fit, rMVPA:::run_global.mvpa_model,
         rMVPA:::.haufe_patterns_from_data),
    function(f) paste(deparse(body(f)), collapse = "\n"),
    character(1)
  )
  expect_false(any(grepl("cov\\(", bodies)))

  imp <- model_importance(fit, X)
  expect_length(imp, P)
  expect_equal(imp, haufe_importance(extract_weights(fit), cov(X))$importance, tolerance = 1e-8)
})

test_that("run_global averaged activation patterns match the dense reference", {
  skip_if_not_installed("sda")
  set.seed(16)
  ds <- gen_sample_dataset(c(4, 4, 4), 40, nlevels = 2, blocks = 4)
  cval <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(load_model("sda_notune"), ds$dataset, ds$design,
                      "classification", crossval = cval)
  res <- run_global(mspec)
  X <- get_feature_matrix(ds$dataset)
  expect_false(is.null(res$activation_patterns))
  ref <- haufe_importance(res$raw_weights, cov(X))$A
  expect_equal(unname(res$activation_patterns), unname(ref), tolerance = 1e-8)
})

test_that("run_global has an informative default for non-mvpa_model specs", {
  fake <- structure(list(), class = c("some_plugin_model", "model_spec"))
  expect_error(run_global(fake), "no global .* method for model class 'some_plugin_model")
})
