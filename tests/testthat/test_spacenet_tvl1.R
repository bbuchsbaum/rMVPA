test_that("spacenet_tvl1 model is registered", {
  mod <- load_model("spacenet_tvl1")
  expect_true(is.list(mod))
  expect_true(is.function(mod$fit))
  expect_true(is.function(mod$predict))
  expect_true(is.function(mod$prob))
})

test_that("spacenet_tvl1 run_global works on binary volumetric data", {
  set.seed(2026)
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("spacenet_tvl1")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)

  res <- run_global(mspec, return_fits = TRUE)

  expect_s3_class(res, "global_mvpa_result")
  expect_true(nrow(res$performance_table) > 0)

  P <- sum(ds$dataset$mask > 0)
  expect_equal(length(res$importance_vector), P)
  expect_true(all(is.finite(res$importance_vector)))
  expect_true(inherits(res$importance_map, "NeuroVol"))

  non_null_fits <- Filter(Negate(is.null), res$fold_fits)
  expect_true(length(non_null_fits) > 0)
  expect_true(inherits(non_null_fits[[1]]$fit, "spacenet_fit"))
})

test_that("spacenet_tvl1 exposes extract_weights and model_importance", {
  set.seed(2027)
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20, nlevels = 2, blocks = 2)
  X <- get_feature_matrix(ds$dataset)
  y <- y_train(ds$design)
  vox <- which(ds$dataset$mask > 0)

  mod <- load_model("spacenet_tvl1")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification",
                      crossval = blocked_cross_validation(ds$design$block_var))
  mfit <- train_model(mspec, as.data.frame(X), y, indices = vox)

  W <- extract_weights(mfit$fit)
  expect_true(is.matrix(W))
  expect_equal(nrow(W), ncol(X))
  expect_equal(ncol(W), 1)

  X_train_masked <- X[, mfit$feature_mask, drop = FALSE]
  imp <- model_importance(mfit$fit, X_train = X_train_masked)
  expect_true(is.numeric(imp))
  expect_equal(length(imp), ncol(X_train_masked))
  expect_true(all(is.finite(imp)))
})

test_that("spacenet PDHG collapses coefficients under extreme l1 penalty", {
  set.seed(2028)
  X <- matrix(rnorm(40 * 6), nrow = 40, ncol = 6)
  y <- rnorm(40)

  fit <- rMVPA:::.spacenet_pdhg_solver(
    X = X,
    y = y,
    alpha = 1e3,
    l1_ratio = 1.0,
    edges = matrix(integer(0), nrow = 0, ncol = 2),
    d_norm2 = 0,
    loss = "mse",
    max_iter = 120,
    tol = 1e-7
  )

  expect_lt(max(abs(fit$w)), 1e-6)
})

test_that("spacenet PDHG is permutation equivariant without TV term", {
  set.seed(2029)
  X <- matrix(rnorm(60 * 5), nrow = 60, ncol = 5)
  y <- rnorm(60)
  perm <- sample.int(ncol(X))

  fit_ref <- rMVPA:::.spacenet_pdhg_solver(
    X = X,
    y = y,
    alpha = 0.05,
    l1_ratio = 0.7,
    edges = matrix(integer(0), nrow = 0, ncol = 2),
    d_norm2 = 0,
    loss = "mse",
    max_iter = 200,
    tol = 1e-6
  )
  fit_perm <- rMVPA:::.spacenet_pdhg_solver(
    X = X[, perm, drop = FALSE],
    y = y,
    alpha = 0.05,
    l1_ratio = 0.7,
    edges = matrix(integer(0), nrow = 0, ncol = 2),
    d_norm2 = 0,
    loss = "mse",
    max_iter = 200,
    tol = 1e-6
  )

  w_perm_back <- numeric(length(perm))
  w_perm_back[perm] <- fit_perm$w
  expect_equal(w_perm_back, fit_ref$w, tolerance = 1e-4)
})

test_that("spacenet core math helpers satisfy contracts", {
  x <- c(-3, -0.2, 0, 0.2, 3)
  thresh <- 0.75
  st <- rMVPA:::.spacenet_soft_threshold(x, thresh)

  expect_true(all(sign(st) == c(-1, 0, 0, 0, 1)))
  expect_equal(abs(st), pmax(abs(x) - thresh, 0))
  expect_equal(rMVPA:::.spacenet_soft_threshold(-x, thresh), -st)

  z <- c(-1200, -2, 0, 2, 1200)
  sig <- rMVPA:::.spacenet_sigmoid(z)
  expect_true(all(sig >= 0 & sig <= 1))
  expect_true(all(is.finite(sig)))
  expect_equal(sig + rMVPA:::.spacenet_sigmoid(-z), rep(1, length(sig)), tolerance = 1e-12)

  eta <- matrix(c(1, 2, 3, 4, -1, -2, 10, -10), nrow = 2L, byrow = TRUE)
  sm <- rMVPA:::.spacenet_softmax(eta)
  expect_equal(rowSums(sm), rep(1, nrow(sm)), tolerance = 1e-12)
  expect_equal(sm, rMVPA:::.spacenet_softmax(eta + 7), tolerance = 1e-12)
})

test_that("spacenet y01 supports factor and numeric encoding contracts", {
  y <- c(-2, -1, 1, 2)
  expect_equal(rMVPA:::.spacenet_y01(y), c(0, 0, 1, 1))

  y_num_const <- c(3, 3, 3)
  expect_equal(rMVPA:::.spacenet_y01(y_num_const), rep(1, 3))

  y_fac <- factor(c("neg", "pos", "neg", "pos"), levels = c("neg", "pos"))
  expect_equal(rMVPA:::.spacenet_y01(y_fac), c(0, 1, 0, 1))
  expect_equal(rMVPA:::.spacenet_y01(y_fac, positive = "neg"), c(1, 0, 1, 0))

  y_single <- factor(rep("a", 4))
  expect_true(all(rMVPA:::.spacenet_y01(y_single) == 0))
})

test_that("spacenet alpha grid is monotone and numerically stable", {
  set.seed(310)
  X <- matrix(rnorm(160), nrow = 20, ncol = 8)
  y <- rbinom(20, size = 1, prob = 0.5)

  alpha_grid <- rMVPA:::.spacenet_make_alpha_grid(
    X = X,
    y = y,
    l1_ratio = 0.7,
    n_alphas = 6L,
    alpha_min_ratio = 1e-3,
    is_classif = TRUE
  )

  expect_length(alpha_grid, 6L)
  expect_true(all(is.finite(alpha_grid)))
  expect_true(all(diff(alpha_grid) < 0))
  expect_true(alpha_grid[1L] >= alpha_grid[length(alpha_grid)])

  deg_grid <- rMVPA:::.spacenet_make_alpha_grid(
    X = matrix(1, nrow = 10L, ncol = 4L),
    y = rnorm(10),
    l1_ratio = 0.5,
    n_alphas = 4L,
    alpha_min_ratio = 1e-3,
    is_classif = FALSE
  )

  expect_true(all(is.finite(deg_grid)))
  expect_true(all(diff(deg_grid) < 0))
  expect_true(deg_grid[1L] >= deg_grid[length(deg_grid)])
})

test_that("spacenet PDHG matches the one-dimensional MSE oracle", {
  set.seed(401)
  n <- 300L
  x <- rnorm(n)
  y <- 2.25 * x + rnorm(n, sd = 0.02)

  X <- matrix(x, ncol = 1L)
  yc <- y - mean(y)
  Xc <- sweep(X, 2L, colMeans(X), "-")
  oracle <- sum(Xc * yc) / max(sum(Xc^2), .Machine$double.eps)
  lspec <- crossprod(Xc, Xc) / max(n, 1L)

  fit <- rMVPA:::.spacenet_pdhg_solver(
    X = Xc,
    y = yc,
    alpha = 0,
    l1_ratio = 1,
    edges = matrix(integer(0), nrow = 0L, ncol = 2L),
    d_norm2 = 0L,
    loss = "mse",
    max_iter = 2000L,
    tol = 1e-8,
    data_lipschitz = lspec
  )

  expect_true(all(is.finite(fit$w)))
  expect_equal(fit$w, oracle, tolerance = 1e-3)
  expect_lte(mean((as.vector(Xc * oracle) - as.vector(Xc * fit$w))^2), 1e-6)
})

test_that("spacenet PDHG preserves permutation invariance with TV coupling", {
  set.seed(402)
  n <- 180L
  p <- 8L
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  y <- rnorm(n)
  edges <- cbind(seq_len(p - 1L), seq_len(p - 1L) + 1L)
  d_norm2 <- 2L * max(tabulate(c(edges), nbins = p), 1L)

  fit_ref <- rMVPA:::.spacenet_pdhg_solver(
    X = X,
    y = y,
    alpha = 0.05,
    l1_ratio = 0.6,
    edges = edges,
    d_norm2 = d_norm2,
    loss = "mse",
    max_iter = 400L,
    tol = 1e-6
  )

  perm <- sample.int(p)
  e1 <- match(edges[, 1L], perm)
  e2 <- match(edges[, 2L], perm)
  keep <- e1 > e2
  e2[keep] <- e1[keep] + e2[keep] - (e1[keep] <- e2[keep])
  edges_perm <- cbind(e1, e2)

  fit_perm <- rMVPA:::.spacenet_pdhg_solver(
    X = X[, perm, drop = FALSE],
    y = y,
    alpha = 0.05,
    l1_ratio = 0.6,
    edges = edges_perm,
    d_norm2 = d_norm2,
    loss = "mse",
    max_iter = 400L,
    tol = 1e-6
  )

  w_perm_back <- numeric(p)
  w_perm_back[perm] <- fit_perm$w
  expect_equal(w_perm_back, fit_ref$w, tolerance = 1e-2)
})

test_that("spacenet CV path short-circuits on small sample sizes", {
  cv_path <- rMVPA:::.spacenet_cv_path(
    X = matrix(rnorm(10), nrow = 5L, ncol = 2L),
    y = rnorm(5),
    feature_ids = 1:2,
    spatial_mask = NULL,
    l1_ratio = 0.5,
    n_alphas = 4L,
    alpha_min_ratio = 1e-3,
    screening_percentile = 30,
    max_iter = 120,
    tol = 1e-4,
    inner_folds = 4L,
    is_classif = FALSE
  )

  expect_length(cv_path$alpha_grid, 4L)
  expect_true(all(is.na(cv_path$cv_loss)))
  expect_equal(cv_path$best_alpha, cv_path$alpha_grid[1L])
})

test_that("spacenet PDHG stays finite on rank-deficient designs", {
  X <- cbind(1, 1 + 1e-8 * (1:24), rep(0, 24))
  y <- rnorm(24)

  fit <- rMVPA:::.spacenet_pdhg_solver(
    X = X,
    y = y,
    alpha = 0.03,
    l1_ratio = 0.75,
    edges = matrix(integer(0), nrow = 0L, ncol = 2L),
    d_norm2 = 0L,
    loss = "mse",
    max_iter = 600L,
    tol = 1e-6
  )

  expect_length(fit$w, ncol(X))
  expect_true(all(is.finite(fit$w)))
  expect_true(all(is.finite(fit$rel_history)))
})

test_that("spacenet_tvl1 fit enforces class cardinality constraints", {
  mod <- load_model("spacenet_tvl1")
  x <- matrix(rnorm(20), nrow = 4L, ncol = 5L)
  y <- factor(rep("A", 4L))
  fit_args <- list(
    l1_ratio = 0.5,
    n_alphas = 6L,
    alpha_min_ratio = 1e-3,
    screening_percentile = 20,
    max_iter = 120L,
    tol = 1e-4,
    inner_folds = 3L
  )

  expect_error(
    mod$fit(x, y, wts = NULL, param = fit_args, lev = levels(y), last = NULL, weights = NULL, classProbs = FALSE),
    "classification requires at least 2 classes"
  )
})

test_that("spacenet_tvl1 supports multibasis global inputs", {
  set.seed(2030)
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 18, nlevels = 2, blocks = 3)
  mb <- mvpa_multibasis_dataset(
    train_data = list(ds$dataset$train_data, ds$dataset$train_data),
    mask = ds$dataset$mask
  )

  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("spacenet_tvl1")
  mspec <- mvpa_model(mod, mb, ds$design, "classification", crossval = cv)
  res <- run_global(mspec)

  P <- sum(ds$dataset$mask > 0) * 2L
  expect_equal(length(res$importance_vector), P)
  expect_true(all(is.finite(res$importance_vector)))
  expect_true(inherits(res$importance_map, "NeuroVol"))
})

test_that("spacenet_tvl1 supports clustered global inputs", {
  set.seed(2031)
  ds <- gen_clustered_sample_dataset(
    D = c(8, 8, 8), nobs = 24, K = 8, nlevels = 2, blocks = 3
  )

  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("spacenet_tvl1")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)
  res <- run_global(mspec)

  P <- ncol(get_feature_matrix(ds$dataset))
  expect_equal(length(res$importance_vector), P)
  expect_true(all(is.finite(res$importance_vector)))
  expect_true(!is.null(res$importance_map))
})

test_that("spacenet_tvl1 supports multiclass global classification", {
  set.seed(2032)
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 30, nlevels = 3, blocks = 3)

  cv <- blocked_cross_validation(ds$design$block_var)
  mod <- load_model("spacenet_tvl1")
  mspec <- mvpa_model(mod, ds$dataset, ds$design, "classification", crossval = cv)
  res <- run_global(mspec, return_fits = TRUE)

  expect_true(nrow(res$performance_table) > 0)
  expect_true(!is.null(res$importance_vector))
  expect_true(!is.null(res$raw_weights))
  expect_true(is.matrix(res$raw_weights))
  expect_equal(ncol(res$raw_weights), 3)

  non_null_fits <- Filter(Negate(is.null), res$fold_fits)
  expect_true(length(non_null_fits) > 0)
  first_fit <- non_null_fits[[1]]$fit
  probs <- mod$prob(first_fit, get_feature_matrix(ds$dataset)[1:5, , drop = FALSE])
  expect_equal(ncol(probs), 3)
  expect_equal(rowSums(probs), rep(1, nrow(probs)), tolerance = 1e-6)
})
