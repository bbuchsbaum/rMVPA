context("banded_ridge_core")

.brc_test_oracle <- function(X, Y, lambdas) {
  x_mu <- colMeans(X)
  x_sd <- apply(X, 2L, stats::sd)
  x_sd[x_sd == 0] <- 1
  Xs <- sweep(sweep(X, 2L, x_mu, "-"), 2L, x_sd, "/")
  y_mu <- colMeans(Y)
  Ys <- sweep(Y, 2L, y_mu, "-")
  beta <- solve(
    crossprod(Xs) + diag(lambdas, nrow = ncol(Xs), ncol = ncol(Xs)),
    crossprod(Xs, Ys)
  )
  coefficients <- sweep(beta, 1L, x_sd, "/")
  intercept <- y_mu - drop(crossprod(x_mu, coefficients))
  list(
    beta = beta,
    coefficients = coefficients,
    intercept = intercept,
    predict = function(newx) sweep(newx %*% coefficients, 2L, intercept, "+")
  )
}

test_that("direct core agrees with an independent generalized-ridge oracle", {
  set.seed(7101)
  X <- matrix(rnorm(36 * 9), 36, 9)
  Y <- matrix(rnorm(36 * 4), 36, 4,
              dimnames = list(NULL, paste0("voxel_", 1:4)))
  groups <- rep(c("low", "semantic", "audio"), c(2, 4, 3))
  lambda_by_group <- c(low = 0.25, semantic = 3, audio = 11)
  penalty <- unname(lambda_by_group[groups])

  fit <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, groups = groups, lambdas = lambda_by_group
  )
  oracle <- .brc_test_oracle(X, Y, penalty)
  got <- rMVPA:::.banded_ridge_predict_direct(fit, X, groups = groups)

  expect_equal(as.vector(fit$standardized_coefficients), as.vector(oracle$beta),
               tolerance = 1e-11)
  expect_equal(as.vector(fit$coefficients), as.vector(oracle$coefficients),
               tolerance = 1e-11)
  expect_equal(unname(fit$intercept), unname(oracle$intercept),
               tolerance = 1e-11)
  expect_equal(unname(got), unname(oracle$predict(X)), tolerance = 1e-11)
  expect_identical(dim(fit$coefficients), c(9L, 4L))
  expect_identical(colnames(fit$coefficients), colnames(Y))
  expect_identical(names(fit$penalty_by_group), c("low", "semantic", "audio"))
  expect_identical(fit$solver$backend, "direct")
  expect_true(all(is.finite(got)))
})

test_that("seeded direct-core property matrix agrees with the independent oracle", {
  set.seed(7102)
  shapes <- list(
    c(n = 8L, p = 1L, v = 1L),
    c(n = 9L, p = 4L, v = 2L),
    c(n = 12L, p = 11L, v = 3L),
    c(n = 7L, p = 13L, v = 4L),
    c(n = 25L, p = 8L, v = 5L)
  )

  for (shape in shapes) {
    for (replicate_id in 1:4) {
      n <- shape[["n"]]
      p <- shape[["p"]]
      v <- shape[["v"]]
      X <- matrix(rnorm(n * p), n, p)
      if (p >= 3L && replicate_id == 4L) {
        X[, p] <- X[, 1] + 1e-5 * rnorm(n)
      }
      Y <- matrix(rnorm(n * v), n, v)
      groups <- rep(c("a", "b", "c"), length.out = p)
      group_names <- unique(groups)
      lambda_by_group <- setNames(
        exp(seq(log(0.05), log(5), length.out = length(group_names))),
        group_names
      )
      penalty <- unname(lambda_by_group[groups])

      fit <- rMVPA:::.banded_ridge_fit_direct(
        X, Y, groups = groups, lambdas = lambda_by_group
      )
      oracle <- .brc_test_oracle(X, Y, penalty)
      got <- rMVPA:::.banded_ridge_predict_direct(fit, X, groups = groups)

      condition_scale <- min(1e-7, 1e-11 * max(1, kappa(
        crossprod(scale(X)) + diag(penalty, nrow = p, ncol = p)
      )))
      expect_equal(unname(got), unname(oracle$predict(X)),
                   tolerance = condition_scale,
                   info = paste(shape, collapse = "x"))
      expect_true(all(is.finite(fit$coefficients)))
    }
  }
})

test_that("equal band penalties reduce to ordinary ridge", {
  set.seed(7103)
  X <- matrix(rnorm(30 * 7), 30, 7)
  Y <- matrix(rnorm(30 * 3), 30, 3)
  groups <- rep(c("vision", "language"), c(3, 4))

  banded <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, groups = groups, lambdas = c(vision = 2.5, language = 2.5)
  )
  ordinary <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, lambdas = c(all = 2.5)
  )

  expect_equal(
    unname(rMVPA:::.banded_ridge_predict_direct(banded, X, groups = groups)),
    unname(rMVPA:::.banded_ridge_predict_direct(ordinary, X)),
    tolerance = 1e-12
  )
})

test_that("feature, band, and response permutations preserve fitted values", {
  set.seed(7104)
  X <- matrix(rnorm(24 * 8), 24, 8)
  Y <- matrix(rnorm(24 * 4), 24, 4)
  groups <- rep(c("a", "b", "c"), c(2, 3, 3))
  lambdas <- c(a = 0.4, b = 2, c = 8)
  base_fit <- rMVPA:::.banded_ridge_fit_direct(X, Y, groups = groups, lambdas = lambdas)
  base_pred <- rMVPA:::.banded_ridge_predict_direct(base_fit, X, groups = groups)

  column_order <- c(7, 3, 1, 8, 5, 2, 6, 4)
  perm_fit <- rMVPA:::.banded_ridge_fit_direct(
    X[, column_order], Y,
    groups = groups[column_order], lambdas = lambdas
  )
  perm_pred <- rMVPA:::.banded_ridge_predict_direct(
    perm_fit, X[, column_order], groups = groups[column_order]
  )
  expect_equal(unname(perm_pred), unname(base_pred), tolerance = 1e-11)

  response_order <- c(4, 2, 1, 3)
  response_fit <- rMVPA:::.banded_ridge_fit_direct(
    X, Y[, response_order], groups = groups, lambdas = lambdas
  )
  response_pred <- rMVPA:::.banded_ridge_predict_direct(response_fit, X, groups = groups)
  expect_equal(unname(response_pred[, order(response_order)]), unname(base_pred),
               tolerance = 1e-12)

  split_pred <- do.call(cbind, lapply(seq_len(ncol(Y)), function(j) {
    fit_j <- rMVPA:::.banded_ridge_fit_direct(
      X, Y[, j], groups = groups, lambdas = lambdas
    )
    rMVPA:::.banded_ridge_predict_direct(fit_j, X, groups = groups)
  }))
  expect_equal(unname(split_pred), unname(base_pred), tolerance = 1e-12)
})

test_that("named-list bands are reordered safely at prediction", {
  set.seed(7105)
  X <- list(
    low = matrix(rnorm(20 * 3), 20, 3,
                 dimnames = list(NULL, paste0("l", 1:3))),
    semantic = matrix(rnorm(20 * 2), 20, 2,
                      dimnames = list(NULL, paste0("s", 1:2)))
  )
  Y <- matrix(rnorm(20 * 2), 20, 2)
  fit <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, lambdas = c(low = 0.5, semantic = 4)
  )
  expected <- rMVPA:::.banded_ridge_predict_direct(fit, X)

  reordered <- rMVPA:::.banded_ridge_predict_direct(fit, X[c("semantic", "low")])
  expect_equal(reordered, expected, tolerance = 0)
  expect_error(
    rMVPA:::.banded_ridge_predict_direct(fit, list(low = X$low)),
    "band names do not match"
  )
  bad_names <- X
  colnames(bad_names$semantic) <- c("other1", "other2")
  expect_error(rMVPA:::.banded_ridge_predict_direct(fit, bad_names),
               "column names/order")
})

test_that("raw feature rescaling is neutral under fold-local standardization", {
  set.seed(7106)
  X <- matrix(rnorm(28 * 6), 28, 6)
  Y <- matrix(rnorm(28 * 3), 28, 3)
  groups <- rep(c("a", "b"), each = 3)
  lambdas <- c(a = 0.7, b = 3)
  base_fit <- rMVPA:::.banded_ridge_fit_direct(X, Y, groups = groups, lambdas = lambdas)
  base_pred <- rMVPA:::.banded_ridge_predict_direct(base_fit, X, groups = groups)

  multipliers <- c(1e-4, 3, 100, 0.2, 7, 1e3)
  X_scaled <- sweep(X, 2, multipliers, "*")
  scaled_fit <- rMVPA:::.banded_ridge_fit_direct(
    X_scaled, Y, groups = groups, lambdas = lambdas
  )
  scaled_pred <- rMVPA:::.banded_ridge_predict_direct(
    scaled_fit, X_scaled, groups = groups
  )
  expect_equal(unname(scaled_pred), unname(base_pred), tolerance = 1e-10)
})

test_that("rank deficiency uses an explicit minimum-norm fallback", {
  set.seed(7107)
  z <- matrix(rnorm(18 * 3), 18, 3)
  X <- cbind(z, z[, 1], z[, 2], 1)
  Y <- matrix(rnorm(18 * 2), 18, 2)
  groups <- rep("all", ncol(X))
  fit <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, groups = groups, lambdas = c(all = 0)
  )

  Xs <- scale(X)
  Xs[, apply(X, 2, stats::sd) == 0] <- 0
  Ys <- scale(Y, center = TRUE, scale = FALSE)
  beta_ref <- MASS::ginv(crossprod(Xs)) %*% crossprod(Xs, Ys)
  pred_ref <- sweep(Xs %*% beta_ref, 2, colMeans(Y), "+")
  pred <- rMVPA:::.banded_ridge_predict_direct(fit, X, groups = groups)

  expect_identical(fit$solver$linear_solver, "direct_symmetric_pseudoinverse")
  expect_lt(fit$solver$rank, ncol(X))
  expect_equal(unname(pred), unname(pred_ref), tolerance = 1e-8)
  expect_true(all(is.finite(pred)))
})

test_that("theta boundary excludes bands without non-finite predictions", {
  set.seed(7108)
  X <- list(
    signal = matrix(rnorm(22 * 3), 22, 3),
    excluded = matrix(rnorm(22 * 4), 22, 4)
  )
  Y <- matrix(rnorm(22 * 2), 22, 2)
  fit <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, alpha = 2, theta = c(signal = 1, excluded = 0)
  )
  signal_only <- rMVPA:::.banded_ridge_fit_direct(
    X$signal, Y, lambdas = c(all = 2)
  )
  pred <- rMVPA:::.banded_ridge_predict_direct(fit, X)
  pred_signal <- rMVPA:::.banded_ridge_predict_direct(signal_only, X$signal)

  expect_identical(fit$excluded_groups, "excluded")
  expect_true(is.infinite(fit$penalty_by_group[["excluded"]]))
  expect_true(all(fit$coefficients[grepl("^excluded::", rownames(fit$coefficients)), ] == 0))
  expect_equal(unname(pred), unname(pred_signal), tolerance = 1e-12)
  expect_true(all(is.finite(pred)))

  near_boundary <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, alpha = 2, theta = c(signal = 1 - 1e-12, excluded = 1e-12)
  )
  expect_true(all(is.finite(rMVPA:::.banded_ridge_predict_direct(near_boundary, X))))

  extreme <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, lambdas = c(signal = 0, excluded = 1e12)
  )
  expect_true(all(is.finite(rMVPA:::.banded_ridge_predict_direct(extreme, X))))
})

test_that("constant predictors and responses have explicit finite behavior", {
  set.seed(7109)
  X <- cbind(constant = 5, varying = rnorm(16))
  Y <- cbind(constant_response = 3, varying_response = rnorm(16))
  fit <- rMVPA:::.banded_ridge_fit_direct(X, Y, lambdas = c(all = 1))
  pred <- rMVPA:::.banded_ridge_predict_direct(fit, X)
  score <- rMVPA:::.banded_ridge_score(Y, pred)

  expect_true(fit$constant_x[["all::constant"]])
  expect_true(all(is.finite(pred)))
  expect_equal(pred[, "constant_response"], rep(3, nrow(X)), tolerance = 1e-12)
  expect_true(is.finite(score$mse[[1]]))
  expect_true(is.na(score$correlation[[1]]))
  expect_true(is.na(score$r2[[1]]))
  expect_true(all(is.finite(unlist(score[2, c("mse", "correlation", "r2")]))))
})

test_that("explicit train/test rows keep preprocessing inside training data", {
  set.seed(7110)
  X <- matrix(rnorm(20 * 4), 20, 4)
  Y <- matrix(rnorm(20 * 2), 20, 2)
  train <- 1:14
  test <- 15:20
  shifted <- X
  shifted[test, ] <- shifted[test, ] + 1e6

  base <- rMVPA:::.banded_ridge_fit_predict_direct(
    X, Y, train, test, lambdas = c(all = 2)
  )
  changed <- rMVPA:::.banded_ridge_fit_predict_direct(
    shifted, Y, train, test, lambdas = c(all = 2)
  )
  expect_identical(base$fit$x_center, changed$fit$x_center)
  expect_identical(base$fit$x_scale, changed$fit$x_scale)
  expect_identical(base$fit$coefficients, changed$fit$coefficients)
  expect_false(isTRUE(all.equal(base$predictions, changed$predictions)))
  expect_identical(rownames(base$predictions), as.character(test))
})

test_that("direct core rejects malformed contracts with specific errors", {
  X <- matrix(1:24, 6, 4)
  Y <- matrix(1:12, 6, 2)

  expect_error(rMVPA:::.banded_ridge_fit_direct(data.frame(X), Y, lambdas = c(all = 1)),
               "numeric matrix or a non-empty named list")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, data.frame(Y), lambdas = c(all = 1)),
               "numeric matrix or numeric vector")
  bad_x <- X; bad_x[1, 1] <- NA_real_
  expect_error(rMVPA:::.banded_ridge_fit_direct(bad_x, Y, lambdas = c(all = 1)),
               "only finite")
  bad_y <- Y; bad_y[1, 1] <- Inf
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, bad_y, lambdas = c(all = 1)),
               "only finite")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y[-1, ], lambdas = c(all = 1)),
               "same number of rows")
  expect_error(rMVPA:::.banded_ridge_fit_direct(
    list(a = matrix(numeric(), nrow = 6, ncol = 0)), Y, lambdas = c(a = 1)
  ), "at least one row and one column")
  expect_error(rMVPA:::.banded_ridge_fit_direct(
    list(X[, 1:2, drop = FALSE], X[, 3:4, drop = FALSE]), Y,
    lambdas = c(a = 1, b = 2)
  ), "list must be named")
  expect_error(rMVPA:::.banded_ridge_fit_direct(
    setNames(list(X[, 1:2, drop = FALSE], X[, 3:4, drop = FALSE]), c("a", "a")), Y,
    lambdas = c(a = 1)
  ), "band names must be unique")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, groups = c("a", "b"),
                                                lambdas = c(a = 1, b = 2)),
               "one value per X column")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, groups = c("a", "", "b", "b"),
                                                lambdas = c(a = 1, b = 2)),
               "non-empty")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, groups = c("a", "a", "b", "b"),
                                                lambdas = c(a = 1)),
               "one value per feature band")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, groups = c("a", "a", "b", "b"),
                                                lambdas = c(a = 1, b = -1)),
               "non-negative")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, lambdas = c(all = 1),
                                                alpha = 1, theta = c(all = 1)),
               "exactly one penalty form")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, alpha = 1, theta = NULL),
               "Both `alpha` and `theta`")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, alpha = -1, theta = c(all = 1)),
               "finite non-negative scalar")
  expect_error(rMVPA:::.banded_ridge_fit_direct(X, Y, alpha = 1, theta = c(all = 0.9)),
               "sum to one")
  expect_error(rMVPA:::.banded_ridge_fit_predict_direct(
    X, Y, train = integer(), test = 1:2, lambdas = c(all = 1)
  ), "non-empty")
  expect_error(rMVPA:::.banded_ridge_fit_predict_direct(
    X, Y, train = c(1, 1, 2), test = 3:4, lambdas = c(all = 1)
  ), "duplicated")
  expect_error(rMVPA:::.banded_ridge_fit_predict_direct(
    X, Y, train = 1:3, test = 3:4, lambdas = c(all = 1)
  ), "disjoint")
  expect_error(rMVPA:::.banded_ridge_fit_predict_direct(
    X, Y, train = 1, test = 2:3, lambdas = c(all = 1)
  ), "at least two")
  expect_error(rMVPA:::.banded_ridge_score(Y, Y[, 1, drop = FALSE]),
               "identical dimensions")
})

test_that("direct core matches the existing alpha_recall zero solve", {
  set.seed(7111)
  X <- matrix(rnorm(26 * 7), 26, 7)
  Y <- matrix(rnorm(26 * 3), 26, 3)
  groups <- rep(c("a", "b"), c(3, 4))
  lambdas <- c(a = 0.75, b = 5)
  pen <- unname(lambdas[groups])

  core <- rMVPA:::.banded_ridge_fit_direct(
    X, Y, groups = groups, lambdas = lambdas
  )
  legacy <- rMVPA:::.br_fit_stacked(
    Xe = X,
    Ye = Y,
    Xr_tr = X[1:2, , drop = FALSE],
    Yr_tr = Y[1:2, , drop = FALSE],
    pen = pen,
    alpha_recall = 0,
    w_rec_tr = c(0, 0)
  )
  legacy_pred <- rMVPA:::.br_predict(X, legacy, mode = "stacked")
  core_pred <- rMVPA:::.banded_ridge_predict_direct(core, X, groups = groups)

  expect_equal(unname(core$standardized_coefficients), unname(legacy$beta),
               tolerance = 1e-11)
  expect_equal(unname(core_pred), unname(legacy_pred), tolerance = 1e-11)
})
