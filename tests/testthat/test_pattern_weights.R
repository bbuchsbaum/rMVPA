library(testthat)

# Observation weights in pattern_model.
#
# The contract under test:
#   * uniform weights are exactly no weights,
#   * integer weights are exactly row replication,
#   * a zero weight is exactly a dropped row (so a poisoned zero-weight row
#     cannot influence the fit),
#   * weights reach the estimator through the design's row_weights or the
#     constructor's `weights` argument, without the old "ignored" warning,
#   * the inner tuning loss is weighted the same way the fit is,
#   * reported metrics stay unweighted.

quiet_run <- function(expr) {
  utils::capture.output(res <- suppressMessages(expr))
  res
}

.sim_xy <- function(n = 48, p = 60, K = 3, snr = 1.5, seed = 1) {
  set.seed(seed)
  y <- factor(rep(letters[seq_len(K)], length.out = n))
  A <- matrix(0, p, K - 1)
  A[1:10, 1] <- rnorm(10) * snr
  A[11:20, 2] <- rnorm(10) * snr
  Yc <- scale(stats::model.matrix(~ y - 1), scale = FALSE)
  S <- crossprod(Yc) / n
  eg <- eigen(S, symmetric = TRUE)
  keep <- eg$values > 1e-8
  Wy <- eg$vectors[, keep, drop = FALSE] %*% diag(1 / sqrt(eg$values[keep]))
  X <- (Yc %*% Wy) %*% t(A) + matrix(rnorm(n * p), n, p)
  list(X = X, y = y, blocks = as.integer(cut(seq_len(n), 3, labels = 1:3)))
}

test_that("weight validation happens at construction and in the core", {
  sim <- sim_pattern_data(n = 45, dims = c(5, 5, 3), K = 3, seed = 21)
  n <- 45
  expect_error(pattern_model(sim$dataset, sim$design, weights = rep(1, n - 1)),
               "44 observation weights but 45")
  expect_error(pattern_model(sim$dataset, sim$design, weights = c(rep(1, n - 1), -1)),
               "finite and non-negative")
  expect_error(pattern_model(sim$dataset, sim$design, weights = c(rep(1, n - 1), NA)),
               "finite and non-negative")
  expect_error(pattern_model(sim$dataset, sim$design, weights = rep(0, n)),
               "positive sum")

  w <- c(rep(2, 15), rep(1, 30))
  spec <- pattern_model(sim$dataset, sim$design, weights = w, rank = 1)
  expect_equal(spec$targets_train$row_weights, w)
  expect_output(print(spec), "observation weights: yes")

  # weights are scale-free in the core
  d <- .sim_xy(seed = 2)
  f1 <- rMVPA:::.pattern_fit(d$X, d$y, rank = 2, weights = rep(c(1, 2), 24))
  f2 <- rMVPA:::.pattern_fit(d$X, d$y, rank = 2, weights = rep(c(10, 20), 24))
  expect_equal(f1$A, f2$A, tolerance = 1e-12)
  expect_equal(f1$noise$D, f2$noise$D, tolerance = 1e-12)
})

test_that("uniform weights reproduce the unweighted fit exactly", {
  d <- .sim_xy(seed = 3)
  base <- rMVPA:::.pattern_fit(d$X, d$y, rank = 2)
  unif <- rMVPA:::.pattern_fit(d$X, d$y, rank = 2, weights = rep(3, nrow(d$X)))
  expect_equal(unif$A, base$A, tolerance = 1e-12)
  expect_equal(unif$C, base$C, tolerance = 1e-12)
  expect_equal(unif$noise$D, base$noise$D, tolerance = 1e-12)
  expect_equal(predict(unif, d$X, type = "prob"), predict(base, d$X, type = "prob"),
               tolerance = 1e-12)
  expect_equal(unif$weights, rep(1, nrow(d$X)))
})

test_that("integer weights are exactly row replication", {
  d <- .sim_xy(n = 36, seed = 4)
  reps <- rep(c(2L, 1L, 1L, 3L), length.out = 36)
  idx <- rep(seq_len(36), times = reps)

  # identity noise so the residual-df convention (which sees n differently for
  # 36 weighted rows and 63 replicated rows) cannot enter the comparison
  ctrl <- pattern_control(max_rank = 2, noise = list(type = "identity"))
  fw <- rMVPA:::.pattern_fit(d$X, d$y, rank = 2, control = ctrl, weights = reps)
  fr <- rMVPA:::.pattern_fit(d$X[idx, , drop = FALSE], d$y[idx], rank = 2, control = ctrl)

  # the whitened target basis is only identified up to rotation/sign, so the
  # comparison is on basis-invariant outputs: probabilities, decoded targets,
  # and fitted (encoded) measurements
  probe <- d$X[1:12, , drop = FALSE]
  expect_equal(predict(fw, probe, type = "prob"), predict(fr, probe, type = "prob"),
               tolerance = 1e-8)
  expect_equal(predict(fw, probe, type = "decode"), predict(fr, probe, type = "decode"),
               tolerance = 1e-8)
  lv <- factor(letters[1:3], levels = levels(d$y))
  ew <- predict(fw, probe, type = "encode", targets = lv)
  er <- predict(fr, probe, type = "encode", targets = lv)
  expect_equal(as.matrix(ew), as.matrix(er), tolerance = 1e-8)

  # weighted class priors match the replicated empirical priors
  tab <- table(d$y[idx]) / length(idx)
  expect_equal(unname(fw$y_transform$priors), unname(as.numeric(tab)), tolerance = 1e-12)
})

test_that("a zero weight is exactly a dropped row, even a poisoned one", {
  d <- .sim_xy(n = 48, seed = 5)
  X <- d$X
  # poison two rows with huge finite garbage; give them weight zero
  X[c(3, 40), ] <- 1e6 * matrix(rnorm(2 * ncol(X)), 2)
  w <- rep(1, 48); w[c(3, 40)] <- 0

  fw <- rMVPA:::.pattern_fit(X, d$y, rank = 2, weights = w)
  fd <- rMVPA:::.pattern_fit(X[-c(3, 40), , drop = FALSE], d$y[-c(3, 40)], rank = 2)

  expect_equal(fw$A, fd$A, tolerance = 1e-12)
  expect_equal(fw$C, fd$C, tolerance = 1e-12)
  expect_equal(fw$noise$D, fd$noise$D, tolerance = 1e-12)
  expect_equal(fw$x_transform$mu, fd$x_transform$mu, tolerance = 1e-12)
  expect_equal(predict(fw, d$X, type = "prob"), predict(fd, d$X, type = "prob"),
               tolerance = 1e-12)
  expect_equal(fw$n_train, 46L)

  # a class carried only by zero-weight rows vanishes from the fit, exactly as
  # if its rows were absent
  w2 <- rep(1, 48); w2[d$y == "c"] <- 0
  f2 <- rMVPA:::.pattern_fit(d$X, d$y, rank = 1, weights = w2, cap_rank = TRUE)
  expect_equal(f2$y_transform$levels, c("a", "b"))

  # fewer than three positively weighted rows is an error
  w3 <- rep(0, 48); w3[1:2] <- 1
  expect_error(rMVPA:::.pattern_fit(d$X, d$y, rank = 1, weights = w3),
               "at least three")
})

test_that("the held-out tuning loss is weighted like the fit", {
  d <- .sim_xy(seed = 6)
  fit <- rMVPA:::.pattern_fit(d$X[1:36, ], d$y[1:36], rank = 2)
  Xte <- d$X[37:48, ]; yte <- d$y[37:48]
  w <- runif(12, 0.2, 2)

  P <- predict(fit, Xte, type = "prob")
  ll <- -log(pmax(P[cbind(seq_len(12), as.integer(yte))], 1e-12))
  expect_equal(rMVPA:::.pattern_loss(fit, Xte, yte, weights = w),
               stats::weighted.mean(ll, w), tolerance = 1e-12)
  # zero-weight assessment rows are excluded; all-zero weights yield NA
  w0 <- w; w0[1:6] <- 0
  expect_equal(rMVPA:::.pattern_loss(fit, Xte, yte, weights = w0),
               stats::weighted.mean(ll[7:12], w0[7:12]), tolerance = 1e-12)
  expect_true(is.na(rMVPA:::.pattern_loss(fit, Xte, yte, weights = rep(0, 12))))
})

test_that("design row_weights flow into fitting without the old warning", {
  sim <- sim_pattern_data(n = 90, dims = c(6, 6, 3), K = NULL, q = 3, snr = 2, seed = 7)
  targets <- sim$targets                       # 90 x 3 continuous targets
  block_var <- sim$design$block_var

  # corrupt block 1: permute its target rows so its X-target pairing is noise
  set.seed(8)
  b1 <- which(block_var == 1)
  targets_bad <- targets
  targets_bad[b1, ] <- targets_bad[sample(b1), ]

  w <- rep(1, 90); w[b1] <- 0
  fs_w <- feature_sets(list(f = targets_bad), row_weights = w)
  fs_u <- feature_sets(list(f = targets_bad))
  make_design <- function(fs) feature_sets_design(X_train = fs, block_var_train = block_var)

  # nonuniform weights no longer warn that they are ignored
  expect_no_warning(pattern_model(sim$dataset, make_design(fs_w), rank = 1))

  spec_w <- pattern_model(sim$dataset, make_design(fs_w), rank = 1)
  spec_u <- pattern_model(sim$dataset, make_design(fs_u), rank = 1)
  res_w <- quiet_run(run_global(spec_w, preflight = "off"))
  res_u <- quiet_run(run_global(spec_u, preflight = "off"))
  expect_true(all(is.finite(unlist(res_w$performance_table))))

  # score only the clean blocks (block 1's shuffled targets are unpredictable
  # for both runs); training that ignores the corrupted block must decode the
  # clean blocks better than training that swallows it
  clean <- function(res) {
    led <- res$ledger
    rows <- which(block_var[led$observation] != 1)
    Y <- as.matrix(led$truth)[rows, , drop = FALSE]
    Yhat <- as.matrix(led$prediction)[rows, , drop = FALSE]
    B <- led$baseline[rows, , drop = FALSE]
    1 - sum((Y - Yhat)^2) / sum((Y - B)^2)
  }
  expect_gt(clean(res_w), clean(res_u))

  # the fold fit actually dropped the zero-weight training rows
  spec_fit <- pattern_model(sim$dataset, make_design(fs_w), rank = 1, return_fits = TRUE)
  res_fit <- quiet_run(run_global(spec_fit, preflight = "off"))
  ff <- res_fit$fold_fits[[which(vapply(res_fit$fold_fits, function(f) f$n_train, 1L) ==
                                   min(vapply(res_fit$fold_fits, function(f) f$n_train, 1L)))[1]]]
  expect_lt(ff$n_train, 60L)                   # a 2-block training set minus block-1 rows
})

test_that("categorical weights work through fit_roi and keep metrics unweighted", {
  sim <- sim_pattern_data(n = 45, dims = c(5, 5, 3), K = 3, seed = 9)
  w <- rep(c(1, 1, 2), 15)
  spec <- pattern_model(sim$dataset, sim$design, rank = 1, weights = w)
  rd <- mock_roi_data(train_data = sim$X[, 1:30], indices = 1:30)
  res <- fit_roi(spec, rd, mock_context(design = spec$design, cv_spec = spec$crossval, id = 1L))
  expect_false(res$error)
  expect_true(is.finite(res$metrics[["Accuracy"]]))

  # uniform weights take the unweighted path end to end: identical metrics
  spec_u <- pattern_model(sim$dataset, sim$design, rank = 1, weights = rep(2, 45))
  spec_0 <- pattern_model(sim$dataset, sim$design, rank = 1)
  res_u <- fit_roi(spec_u, rd, mock_context(design = spec_u$design, cv_spec = spec_u$crossval, id = 1L))
  res_0 <- fit_roi(spec_0, rd, mock_context(design = spec_0$design, cv_spec = spec_0$crossval, id = 1L))
  expect_equal(res_u$metrics, res_0$metrics, tolerance = 1e-12)
})
