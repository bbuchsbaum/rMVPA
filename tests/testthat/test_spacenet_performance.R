test_that("spacenet PDHG has bounded iteration growth on structured data", {
  set.seed(901)
  n <- 200L
  p <- 60L
  x_raw <- qr.Q(qr(matrix(rnorm(n * p), nrow = n, ncol = p)))

  w_true <- numeric(p)
  w_true[1:6] <- c(1.2, -0.9, 0.7, -0.6, 0.4, 0.3)
  y <- as.vector(x_raw %*% w_true + rnorm(n, sd = 0.05))
  y <- y - mean(y)

  fit <- rMVPA:::.spacenet_pdhg_solver(
    X = x_raw,
    y = y,
    alpha = 0.005,
    l1_ratio = 0.8,
    edges = matrix(integer(0), nrow = 0L, ncol = 2L),
    d_norm2 = 0L,
    loss = "mse",
    max_iter = 600L,
    tol = 1e-4,
    data_lipschitz = rMVPA:::.spacenet_spectral_norm_squared(x_raw) / n
  )

  expect_true(is.finite(fit$rel_history[length(fit$rel_history)]))
  expect_true(all(is.finite(fit$rel_history)))
  expect_lte(length(fit$rel_history), 600L)
  expect_true(length(fit$rel_history) <= 400L)
  expect_lt(fit$rel_history[length(fit$rel_history)], 5e-3)
  expect_true(fit$rel_history[1L] > fit$rel_history[length(fit$rel_history)] * 5)
})
