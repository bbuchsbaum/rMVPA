context("banded_ridge_solvers")

.brs_test_predictions <- function(path, X, groups) {
  lapply(path$fits, function(fit) {
    rMVPA:::.banded_ridge_predict_optimized(fit, X, groups = groups)
  })
}

.brs_test_compare <- function(X, Y, groups, theta, alphas, newx = X) {
  svd_path <- rMVPA:::.banded_ridge_fit_svd_path(
    X, Y, alphas, theta, groups = groups, recover_dual = TRUE
  )
  dual_path <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, alphas, theta, groups = groups, recover_primal = TRUE
  )
  direct <- lapply(alphas, function(alpha) {
    rMVPA:::.banded_ridge_fit_direct(
      X, Y, groups = groups, alpha = alpha, theta = theta
    )
  })
  list(
    direct = direct,
    svd = svd_path,
    dual = dual_path,
    direct_predictions = lapply(direct, function(fit) {
      rMVPA:::.banded_ridge_predict_direct(fit, newx, groups = groups)
    }),
    svd_predictions = .brs_test_predictions(svd_path, newx, groups),
    dual_predictions = .brs_test_predictions(dual_path, newx, groups)
  )
}

test_that("SVD and dual backends agree with direct primal for p<n, p=n, and p>n", {
  shapes <- list(c(n = 30L, p = 8L, v = 1L),
                 c(n = 18L, p = 18L, v = 3L),
                 c(n = 12L, p = 31L, v = 5L))
  for (shape_id in seq_along(shapes)) {
    shape <- shapes[[shape_id]]
    set.seed(7900L + shape_id)
    X <- matrix(rnorm(shape[["n"]] * shape[["p"]]), shape[["n"]], shape[["p"]])
    Y <- matrix(rnorm(shape[["n"]] * shape[["v"]]), shape[["n"]], shape[["v"]],
                dimnames = list(NULL, paste0("v", seq_len(shape[["v"]]))))
    groups <- rep(c("low", "semantic", "audio"), length.out = shape[["p"]])
    theta <- c(low = 0.15, semantic = 0.35, audio = 0.5)
    compared <- .brs_test_compare(X, Y, groups, theta, c(0.03, 1.5, 80))

    expect_identical(compared$svd$decomposition_count, 1L)
    expect_identical(compared$dual$decomposition_count, 1L)
    for (aa in seq_along(compared$direct)) {
      scale <- max(1, kappa(crossprod(scale(X)) + diag(
        (c(0.03, 1.5, 80)[[aa]] / theta[groups]), ncol(X), ncol(X)
      )))
      tolerance <- min(1e-7, 5e-12 * scale)
      expect_equal(unname(compared$svd_predictions[[aa]]),
                   unname(compared$direct_predictions[[aa]]), tolerance = tolerance)
      expect_equal(unname(compared$dual_predictions[[aa]]),
                   unname(compared$direct_predictions[[aa]]), tolerance = tolerance)
      expect_equal(unname(compared$svd$fits[[aa]]$standardized_coefficients),
                   unname(compared$direct[[aa]]$standardized_coefficients),
                   tolerance = tolerance)
      expect_equal(unname(compared$dual$fits[[aa]]$standardized_coefficients),
                   unname(compared$direct[[aa]]$standardized_coefficients),
                   tolerance = tolerance)
      expect_equal(compared$svd$fits[[aa]]$intercept,
                   compared$direct[[aa]]$intercept, tolerance = tolerance)
      expect_equal(compared$dual$fits[[aa]]$intercept,
                   compared$direct[[aa]]$intercept, tolerance = tolerance)
    }
  }
})

test_that("seeded differential matrix covers rank, conditioning, scales, and response counts", {
  set.seed(7904)
  shapes <- list(c(n = 8L, p = 1L, v = 1L), c(n = 20L, p = 7L, v = 4L),
                 c(n = 10L, p = 22L, v = 2L), c(n = 16L, p = 16L, v = 6L))
  max_errors <- c(svd_prediction = 0, dual_prediction = 0,
                  svd_weight = 0, dual_weight = 0)
  for (shape_id in seq_along(shapes)) {
    shape <- shapes[[shape_id]]
    for (replicate_id in 1:5) {
      X <- matrix(rnorm(shape[["n"]] * shape[["p"]]), shape[["n"]], shape[["p"]])
      if (shape[["p"]] >= 4L && replicate_id == 2L) X[, shape[["p"]]] <- X[, 1]
      if (shape[["p"]] >= 4L && replicate_id == 3L) {
        X[, shape[["p"]]] <- X[, 1] + 1e-9 * rnorm(shape[["n"]])
      }
      X <- sweep(X, 2L, 10^seq(-6, 6, length.out = shape[["p"]]), "*")
      Y <- matrix(rnorm(shape[["n"]] * shape[["v"]]), shape[["n"]], shape[["v"]])
      groups <- rep(c("a", "b"), length.out = shape[["p"]])
      theta <- if (shape[["p"]] == 1L) c(a = 1) else c(a = 0.23, b = 0.77)
      if (shape[["p"]] == 1L) groups <- "a"
      compared <- .brs_test_compare(X, Y, groups, theta, 0.7)
      direct <- compared$direct_predictions[[1]]
      errors <- c(
        svd_prediction = max(abs(compared$svd_predictions[[1]] - direct)),
        dual_prediction = max(abs(compared$dual_predictions[[1]] - direct)),
        svd_weight = max(abs(compared$svd$fits[[1]]$standardized_coefficients -
                                   compared$direct[[1]]$standardized_coefficients)),
        dual_weight = max(abs(compared$dual$fits[[1]]$standardized_coefficients -
                                    compared$direct[[1]]$standardized_coefficients))
      )
      max_errors <- pmax(max_errors, errors)
      expect_true(all(errors < 2e-7), info = paste(shape_id, replicate_id))
      expect_true(all(is.finite(compared$dual_predictions[[1]])))
    }
  }
  expect_true(all(is.finite(max_errors)))
})

test_that("weighted dual kernel equals an explicit standardized primal and kernel solve", {
  set.seed(7905)
  X <- matrix(rnorm(22 * 9), 22, 9)
  Y <- matrix(rnorm(22 * 3), 22, 3,
              dimnames = list(NULL, c("v1", "v2", "v3")))
  newx <- matrix(rnorm(7 * 9), 7, 9)
  groups <- rep(c("a", "b", "c"), c(2, 3, 4))
  theta <- c(a = 0.1, b = 0.3, c = 0.6)
  alpha <- 2.4
  dual <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, alpha, theta, groups = groups, recover_primal = TRUE
  )$fits[[1]]
  got <- rMVPA:::.banded_ridge_predict_optimized(dual, newx, groups = groups)

  x_center <- colMeans(X)
  x_scale <- apply(X, 2, stats::sd); x_scale[x_scale == 0] <- 1
  Xs <- sweep(sweep(X, 2, x_center, "-"), 2, x_scale, "/")
  Xnew_s <- sweep(sweep(newx, 2, x_center, "-"), 2, x_scale, "/")
  y_center <- colMeans(Y)
  Yc <- sweep(Y, 2, y_center, "-")
  K <- Reduce(`+`, lapply(names(theta), function(g) {
    idx <- groups == g
    theta[[g]] * tcrossprod(Xs[, idx, drop = FALSE])
  }))
  Ktest <- Reduce(`+`, lapply(names(theta), function(g) {
    idx <- groups == g
    theta[[g]] * Xnew_s[, idx, drop = FALSE] %*% t(Xs[, idx, drop = FALSE])
  }))
  weights <- solve(K + diag(alpha, nrow(K)), Yc)
  expected <- sweep(Ktest %*% weights, 2, y_center, "+")
  expect_equal(unname(got), unname(expected), tolerance = 1e-11)

  penalty <- unname(alpha / theta[groups])
  beta <- solve(crossprod(Xs) + diag(penalty, ncol(Xs)), crossprod(Xs, Yc))
  primal_expected <- sweep(Xnew_s %*% beta, 2, y_center, "+")
  expect_equal(unname(got), unname(primal_expected), tolerance = 1e-11)
})

test_that("dual cross-kernel centering is derived only from training rows", {
  set.seed(7906)
  X <- matrix(rnorm(24 * 6), 24, 6)
  Y <- matrix(rnorm(24 * 2), 24, 2,
              dimnames = list(NULL, c("v1", "v2")))
  test_x <- matrix(rnorm(8 * 6), 8, 6)
  groups <- rep(c("a", "b"), each = 3)
  theta <- c(a = 0.4, b = 0.6)
  fit <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, 1, theta, groups = groups
  )$fits[[1]]
  training_state <- list(fit$x_center, fit$x_scale, fit$training_xs,
                         fit$dual_weights)
  base <- rMVPA:::.banded_ridge_predict_optimized(fit, test_x, groups = groups)
  shifted <- sweep(test_x, 2, 10^(3:8), "+")
  changed <- rMVPA:::.banded_ridge_predict_optimized(fit, shifted, groups = groups)
  expect_identical(training_state,
                   list(fit$x_center, fit$x_scale, fit$training_xs,
                        fit$dual_weights))
  expect_false(isTRUE(all.equal(base, changed)))

  Xs_shifted <- sweep(sweep(shifted, 2, fit$x_center, "-"), 2, fit$x_scale, "/")
  explicit_kernel <- Reduce(`+`, lapply(names(theta), function(g) {
    idx <- groups == g
    theta[[g]] * Xs_shifted[, idx, drop = FALSE] %*%
      t(fit$training_xs[, idx, drop = FALSE])
  }))
  explicit <- sweep(explicit_kernel %*% fit$dual_weights, 2, fit$y_center, "+")
  expect_equal(unname(changed), unname(explicit), tolerance = 1e-12)

  leaky_center <- colMeans(rbind(X, shifted))
  leaky_scale <- apply(rbind(X, shifted), 2, stats::sd)
  wrong_train <- sweep(sweep(X, 2, leaky_center, "-"), 2, leaky_scale, "/")
  wrong_test <- sweep(sweep(shifted, 2, leaky_center, "-"), 2, leaky_scale, "/")
  wrong_kernel <- Reduce(`+`, lapply(names(theta), function(g) {
    idx <- groups == g
    theta[[g]] * wrong_test[, idx, drop = FALSE] %*%
      t(wrong_train[, idx, drop = FALSE])
  }))
  wrong <- sweep(wrong_kernel %*% fit$dual_weights, 2, fit$y_center, "+")
  expect_false(isTRUE(all.equal(changed, wrong)))
})

test_that("cache instrumentation reuses only identical fold and theta decompositions", {
  set.seed(7907)
  X <- matrix(rnorm(30 * 8), 30, 8)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  groups <- rep(c("a", "b"), each = 4)
  folds <- rMVPA:::.banded_ridge_make_folds(1:30, 3L, seed = 7908L)
  thetas <- list(c(a = 0.2, b = 0.8), c(a = 0.75, b = 0.25))
  cache <- rMVPA:::.banded_ridge_solver_cache()
  for (ff in seq_along(folds)) {
    for (tt in seq_along(thetas)) {
      rows <- folds[[ff]]$train
      key <- paste0("fold", ff, "-theta", tt)
      first <- rMVPA:::.banded_ridge_fit_svd_path(
        X[rows, ], Y[rows, ], c(0.1, 1, 10), thetas[[tt]], groups = groups,
        cache = cache, cache_key = key
      )
      second <- rMVPA:::.banded_ridge_fit_svd_path(
        X[rows, ], Y[rows, ] * 2, c(0.2, 2), thetas[[tt]], groups = groups,
        cache = cache, cache_key = key
      )
      expect_identical(first$decomposition_count, 1L)
      expect_identical(second$decomposition_count, 0L)
      expect_true(second$cache_hit)
    }
  }
  expect_identical(cache$decomposition_count,
                   as.integer(length(folds) * length(thetas)))
  expect_identical(cache$hit_count,
                   as.integer(length(folds) * length(thetas)))
  expect_error(rMVPA:::.banded_ridge_fit_svd_path(
    X[folds[[1]]$train, ], Y[folds[[1]]$train, ], 1, c(a = 0.5, b = 0.5),
    groups = groups, cache = cache, cache_key = "fold1-theta1"
  ), "cache_key collision")
  expect_error(rMVPA:::.banded_ridge_fit_svd_path(
    X[folds[[1]]$train, ], Y[folds[[1]]$train, ], 1, thetas[[1]],
    groups = groups, solver_tolerance = 1e-8,
    cache = cache, cache_key = "fold1-theta1"
  ), "cache_key collision")

  dual_cache <- rMVPA:::.banded_ridge_solver_cache()
  first_dual <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, c(0.1, 1, 10), thetas[[1]], groups = groups,
    cache = dual_cache, cache_key = "all-theta1"
  )
  second_dual <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, 4, thetas[[1]], groups = groups,
    cache = dual_cache, cache_key = "all-theta1"
  )
  expect_identical(first_dual$decomposition_count, 1L)
  expect_identical(second_dual$decomposition_count, 0L)
  expect_identical(dual_cache$decomposition_count, 1L)
})

test_that("band, response, alpha-batch, and response-batch permutations preserve outputs", {
  set.seed(7909)
  X <- matrix(rnorm(26 * 10), 26, 10)
  Y <- matrix(rnorm(26 * 5), 26, 5,
              dimnames = list(NULL, paste0("v", 1:5)))
  groups <- rep(c("a", "b", "c"), c(2, 3, 5))
  theta <- c(a = 0.15, b = 0.25, c = 0.6)
  alphas <- c(0.02, 0.5, 7, 100)
  base_svd <- rMVPA:::.banded_ridge_fit_svd_path(
    X, Y, alphas, theta, groups = groups,
    alpha_batch_size = 1L, response_batch_size = 1L
  )
  batched_svd <- rMVPA:::.banded_ridge_fit_svd_path(
    X, Y, alphas, theta, groups = groups,
    alpha_batch_size = Inf, response_batch_size = Inf
  )
  base_dual <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, alphas, theta, groups = groups,
    alpha_batch_size = 2L, response_batch_size = 2L
  )
  batched_dual <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, alphas, theta, groups = groups,
    alpha_batch_size = Inf, response_batch_size = Inf
  )
  for (aa in seq_along(alphas)) {
    expect_equal(
      rMVPA:::.banded_ridge_predict_optimized(base_svd$fits[[aa]], X, groups),
      rMVPA:::.banded_ridge_predict_optimized(batched_svd$fits[[aa]], X, groups),
      tolerance = 1e-14
    )
    expect_equal(
      rMVPA:::.banded_ridge_predict_optimized(base_dual$fits[[aa]], X, groups),
      rMVPA:::.banded_ridge_predict_optimized(batched_dual$fits[[aa]], X, groups),
      tolerance = 1e-14
    )
  }

  column_order <- c(10, 3, 1, 8, 5, 2, 6, 4, 9, 7)
  permuted <- rMVPA:::.banded_ridge_fit_dual_path(
    X[, column_order], Y, alphas, theta, groups = groups[column_order]
  )
  for (aa in seq_along(alphas)) {
    base_pred <- rMVPA:::.banded_ridge_predict_optimized(
      base_dual$fits[[aa]], X, groups
    )
    perm_pred <- rMVPA:::.banded_ridge_predict_optimized(
      permuted$fits[[aa]], X[, column_order], groups[column_order]
    )
    expect_equal(unname(base_pred), unname(perm_pred), tolerance = 1e-11)
  }

  response_order <- c(5, 2, 4, 1, 3)
  response_fit <- rMVPA:::.banded_ridge_fit_svd_path(
    X, Y[, response_order], alphas, theta, groups = groups
  )
  for (aa in seq_along(alphas)) {
    expected <- rMVPA:::.banded_ridge_predict_optimized(
      base_svd$fits[[aa]], X, groups
    )
    got <- rMVPA:::.banded_ridge_predict_optimized(
      response_fit$fits[[aa]], X, groups
    )
    expect_equal(unname(got[, order(response_order)]), unname(expected),
                 tolerance = 1e-12)
  }
})

test_that("alpha/theta and explicit generalized penalties name the same estimator", {
  set.seed(7910)
  X <- matrix(rnorm(28 * 7), 28, 7)
  Y <- matrix(rnorm(28 * 3), 28, 3)
  groups <- rep(c("a", "b"), c(2, 5))
  theta <- c(a = 0.3, b = 0.7)
  alpha <- 4
  optimized <- rMVPA:::.banded_ridge_fit_svd_path(
    X, Y, alpha, theta, groups = groups
  )$fits[[1]]
  explicit <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, groups = groups, lambdas = alpha / theta
  )
  expect_equal(
    rMVPA:::.banded_ridge_predict_optimized(optimized, X, groups),
    rMVPA:::.banded_ridge_predict_direct(explicit, X, groups),
    tolerance = 1e-12
  )
  expect_equal(optimized$penalty_by_group, alpha / theta, tolerance = 0)
})

test_that("zero-weight bands, duplicate kernels, singularity, and extreme alpha remain explicit", {
  set.seed(7911)
  base <- matrix(rnorm(8 * 3), 8, 3)
  X <- cbind(base, base, 1, base[, 1])
  Y <- matrix(rnorm(8 * 2), 8, 2,
              dimnames = list(NULL, c("v1", "v2")))
  groups <- c(rep("kept", 3), rep("duplicate", 3), "constant", "kept")
  theta <- c(kept = 1, duplicate = 0, constant = 0)
  alphas <- c(0, 1e-10, 1e10)
  compared <- .brs_test_compare(X, Y, groups, theta, alphas)
  for (aa in seq_along(alphas)) {
    stress_tolerance <- if (alphas[[aa]] == 1e-10) 1e-4 else 1e-7
    expect_true(all(is.finite(compared$svd_predictions[[aa]])))
    expect_true(all(is.finite(compared$dual_predictions[[aa]])))
    expect_equal(unname(compared$svd_predictions[[aa]]),
                 unname(compared$direct_predictions[[aa]]),
                 tolerance = stress_tolerance)
    expect_equal(unname(compared$dual_predictions[[aa]]),
                 unname(compared$direct_predictions[[aa]]),
                 tolerance = stress_tolerance)
    expect_identical(compared$svd$fits[[aa]]$solver$jitter, 0)
    expect_false(compared$svd$fits[[aa]]$solver$estimator_change)
  }
  expect_identical(compared$svd$fits[[1]]$solver$spectral_fallback,
                   "minimum_norm_pseudoinverse")
  expect_identical(compared$dual$fits[[1]]$solver$spectral_fallback,
                   "minimum_norm_pseudoinverse")
  expect_true(all(compared$svd$fits[[2]]$coefficients[
    grepl("^(duplicate|constant)::", rownames(compared$svd$fits[[2]]$coefficients)),
  ] == 0))
})

test_that("dual prediction works without materializing primal coefficients", {
  set.seed(7912)
  X <- matrix(rnorm(12 * 40), 12, 40)
  Y <- matrix(rnorm(12 * 4), 12, 4,
              dimnames = list(NULL, paste0("v", 1:4)))
  groups <- rep(c("a", "b"), each = 20)
  theta <- c(a = 0.4, b = 0.6)
  compact <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, 2, theta, groups = groups, recover_primal = FALSE
  )$fits[[1]]
  full <- rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, 2, theta, groups = groups, recover_primal = TRUE
  )$fits[[1]]
  expect_null(compact$coefficients)
  expect_null(compact$standardized_coefficients)
  expect_null(compact$intercept)
  expect_true(is.matrix(compact$dual_weights))
  expect_equal(
    rMVPA:::.banded_ridge_predict_optimized(compact, X, groups),
    rMVPA:::.banded_ridge_predict_optimized(full, X, groups),
    tolerance = 1e-12
  )
  expect_lt(as.numeric(object.size(compact)), as.numeric(object.size(full)))
  expect_identical(compact$solver$peak_intermediate_dimensions[["primal_batch"]], 0)
})

test_that("solver planning is conservative, overrideable, and memory bounded", {
  small <- rMVPA:::.banded_ridge_solver_plan(
    n = 100, p = 50, v = 20, n_bands = 3, n_alphas = 1
  )
  path <- rMVPA:::.banded_ridge_solver_plan(
    n = 100, p = 100, v = 20, n_bands = 3, n_alphas = 8
  )
  wide <- rMVPA:::.banded_ridge_solver_plan(
    n = 50, p = 500, v = 20, n_bands = 3, n_alphas = 8
  )
  expect_identical(small$solver, "direct")
  expect_match(small$reason, "p <= 64")
  expect_identical(path$solver, "svd_primal")
  expect_identical(wide$solver, "dual_kernel")
  expect_match(wide$reason, "p >= 4n")
  expect_true(wide$fits_memory_contract[[wide$solver]])
  expect_lte(wide$estimated_mib[[wide$solver]], wide$memory_limit_mb)

  override <- rMVPA:::.banded_ridge_solver_plan(
    n = 50, p = 500, v = 20, n_bands = 3, n_alphas = 8,
    solver = "svd_primal", memory_limit_mb = 1024
  )
  expect_identical(override$solver, "svd_primal")
  expect_match(override$reason, "explicit override")
  expect_error(rMVPA:::.banded_ridge_solver_plan(
    n = 500, p = 5000, v = 100, n_bands = 20, n_alphas = 20,
    alpha_batch_size = 10, response_batch_size = 100,
    solver = "dual_kernel", memory_limit_mb = 1
  ), "exceeds")
  expect_error(rMVPA:::.banded_ridge_solver_plan(
    n = 500, p = 5000, v = 100, n_bands = 20, n_alphas = 20,
    alpha_batch_size = 10, response_batch_size = 100,
    memory_limit_mb = 0.001
  ), "No solver fits")
})

test_that("solver contracts reject invalid paths, theta, caches, and prediction layouts", {
  X <- matrix(1:24, 6, 4)
  Y <- matrix(1:12, 6, 2)
  groups <- c("a", "a", "b", "b")
  expect_error(rMVPA:::.banded_ridge_fit_svd_path(
    X, Y, -1, c(a = 0.5, b = 0.5), groups
  ), "non-negative")
  expect_error(rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, 1, c(a = 0.2, b = 0.2), groups
  ), "sum to one")
  expect_error(rMVPA:::.banded_ridge_fit_svd_path(
    X, Y, 1, c(a = 0.5, b = 0.5), groups, response_batch_size = 0
  ), "positive integer or Inf")
  cache <- rMVPA:::.banded_ridge_solver_cache()
  expect_error(rMVPA:::.banded_ridge_fit_dual_path(
    X, Y, 1, c(a = 0.5, b = 0.5), groups, cache = cache
  ), "cache_key")
  fit <- rMVPA:::.banded_ridge_fit_svd_path(
    X, Y, 1, c(a = 0.5, b = 0.5), groups
  )$fits[[1]]
  expect_error(rMVPA:::.banded_ridge_predict_optimized(fit, X[, -1], groups[-1]),
               "group membership/order|different number of columns")
  expect_error(rMVPA:::.banded_ridge_predict_optimized(list(), X, groups),
               "banded_ridge_optimized_fit")
})
