confirmation_fixture <- function(seed = 550, n = 100, p = 5) {
  set.seed(seed)
  Y <- cbind(a = rnorm(n), b = rnorm(n))
  X <- Y %*% matrix(c(1, 0, 0, 1, 1, -1, 0, 0, 0.5, 0.3), 2, p) + matrix(rnorm(n * p), n)
  fit <- .pattern_fit(X, Y, rank = 2, control = pattern_control())
  Ytest <- cbind(a = rnorm(n), b = rnorm(n))
  Xtest <- Ytest %*% matrix(c(1, 0, 0, 1, 1, -1, 0, 0, 0.5, 0.3), 2, p) + matrix(rnorm(n * p), n)
  list(fit = fit, X = X, Y = Y, Xtest = Xtest, Ytest = Ytest,
       discovery = paste0("s1:discovery:", seq_len(n)), ids = paste0("s1:confirm:", seq_len(n)))
}
confirm_fixture <- function(d, ...) pattern_confirm(d$fit, d$Xtest, d$Ytest,
  observation_ids = d$ids, discovery_ids = d$discovery,
  feature_ids = paste0("v", seq_len(ncol(d$X))), preprocessing_id = "raw-BOLD-v1", ...)

test_that("confirmation matches independent lm coefficients and tests", {
  d <- confirmation_fixture(); N <- cbind(drift = seq_len(nrow(d$Xtest)))
  c <- confirm_fixture(d, nuisance = N)
  Tm <- .pattern_targets_apply(d$fit$y_transform, d$Ytest) %*% d$fit$C
  for (v in 1:5) {
    ref <- stats::lm(d$Xtest[, v] ~ Tm + N)
    sm <- summary(ref)
    expect_equal(unname(c$estimate[v, ]), unname(coef(ref)[2:3]), tolerance = 1e-10)
    expect_equal(unname(c$se[v, ]), unname(sm$coefficients[2:3, 2]), tolerance = 1e-10)
    expect_equal(.pattern_confirmation_cov(c, v), unname(vcov(ref)[2:3, 2:3]), tolerance = 1e-10)
    reduced <- lm(d$Xtest[, v] ~ N)
    expect_equal(c$omnibus$F[v], anova(reduced, ref)$F[2], tolerance = 1e-10)
  }
  expect_equal(as.vector(c$p_holm), p.adjust(as.vector(c$p), "holm"))
  expect_equal(c$basis$matrix, pattern_basis(d$fit)$matrix)
})

test_that("CR1 covariance matches a dense cluster sandwich oracle", {
  d <- confirmation_fixture(); blocks <- rep(1:10, each = 10)
  c <- confirm_fixture(d, block_var = blocks, inference = confirmation_plan("block_robust"))
  Tm <- .pattern_targets_apply(d$fit$y_transform, d$Ytest) %*% d$fit$C
  D <- cbind(Tm, 1); bread <- solve(crossprod(D)); beta <- bread %*% crossprod(D, d$Xtest)
  E <- d$Xtest - D %*% beta
  for (v in 1:5) {
    scores <- rowsum(D * E[, v], blocks)
    V <- bread %*% crossprod(scores) %*% bread * (10/9) * (99/97)
    expect_equal(.pattern_confirmation_cov(c, v), V[1:2, 1:2], tolerance = 1e-10)
  }
  expect_equal(c$df, 9)
})

test_that("omnibus and covariance respect frozen coordinate rotations", {
  d <- confirmation_fixture(); c <- confirm_fixture(d)
  Q <- matrix(c(cos(0.6), sin(0.6), -sin(0.6), cos(0.6)), 2)
  rotated <- d; rotated$fit$C <- d$fit$C %*% Q
  cr <- confirm_fixture(rotated)
  expect_equal(unname(cr$estimate), unname(c$estimate %*% Q), tolerance = 1e-10)
  expect_equal(cr$omnibus, c$omnibus, tolerance = 1e-10)
  expect_equal(.pattern_confirmation_cov(cr, 1), crossprod(Q, .pattern_confirmation_cov(c, 1) %*% Q), tolerance = 1e-10)
  expect_false(identical(c$basis$basis_id, cr$basis$basis_id))
})

test_that("provenance, nuisance, and feature contracts fail closed", {
  d <- confirmation_fixture()
  bad <- d; bad$ids[1] <- d$discovery[1]
  expect_error(confirm_fixture(bad), "overlap")
  bad <- d; bad$ids[1] <- bad$ids[2]
  expect_error(confirm_fixture(bad), "unique")
  expect_error(confirm_fixture(d, nuisance = matrix(1, 100, 1)), "rank deficient")
  expect_error(confirm_fixture(d, block_var = rep(1:3, length.out = 100), inference = confirmation_plan("block_robust")), "independent blocks")
  expect_error(confirm_fixture(d, block_var = 1:100), "Independent inference")
  bad <- d; bad$Ytest <- bad$Ytest[, 2:1]
  expect_error(confirm_fixture(bad), "response IDs")
  bad <- d; bad$Xtest[1] <- NA
  expect_error(confirm_fixture(bad), "finite numeric")
  expect_error(confirmation_plan(n_resamples = 2), "at least 99")
  expect_error(confirmation_plan(seed = -1), "nonnegative")
})

test_that("wild bootstrap is reproducible, bounded, and preserves RNG", {
  d <- confirmation_fixture(n = 60); blocks <- rep(1:10, each = 6)
  plan <- confirmation_plan("sign_flip", 99, 314)
  before <- .Random.seed
  c <- confirm_fixture(d, block_var = blocks, inference = plan)
  expect_identical(.Random.seed, before)
  expect_identical(c, confirm_fixture(d, block_var = blocks, inference = plan))
  expect_true(all(c$bootstrap$p_max >= c$p))
  expect_true(all(c$bootstrap$p_omnibus_max >= c$omnibus$p))
  expect_true(all(c$p >= 0.01 & c$p <= 1))
  expect_equal(c$estimate, confirm_fixture(d)$estimate)
})

test_that("component heads are calibrated without confirmation leakage", {
  d <- confirmation_fixture(); c <- confirm_fixture(d)
  cal <- list(dataset = d$X, design = d$Y, observation_ids = d$discovery)
  out <- pattern_component_tests(d$fit, c, d$Xtest, d$Ytest, cal)
  Z <- predict(d$fit, d$X, type = "scores")
  Ztest <- predict(d$fit, d$Xtest, type = "scores")
  expected <- lm(d$Y ~ Z)
  expect_equal(unname(out$heads$full), unname(coef(expected)), tolerance = 1e-10)
  expect_equal(out$losses$full, rowMeans((d$Ytest - cbind(1, Ztest) %*% coef(expected))^2))
  for (k in 1:2) {
    delta <- out$losses$reduced[, k] - out$losses$full
    ref <- t.test(delta)
    expect_equal(out$incremental$p[k], ref$p.value, tolerance = 1e-10)
    Tm <- .pattern_targets_apply(d$fit$y_transform, d$Ytest) %*% d$fit$C
    expect_equal(out$association$p[k], cor.test(Ztest[, k], Tm[, k])$p.value, tolerance = 1e-10)
  }
  expect_error(pattern_component_tests(d$fit, c, d$Xtest + 1, d$Ytest, cal), "must match")
  bad <- cal; bad$observation_ids[1] <- d$ids[1]
  expect_error(pattern_component_tests(d$fit, c, d$Xtest, d$Ytest, bad), "overlap")
  expect_null(pattern_component_tests(d$fit, c, d$Xtest, d$Ytest)$incremental)
})

test_that("independent null p values calibrate and signal has power across seeds", {
  d <- confirmation_fixture(); Tm <- .pattern_targets_apply(d$fit$y_transform, d$Ytest) %*% d$fit$C
  s <- .pattern_lm_setup(cbind(Tm, 1), 2)
  set.seed(555)
  # Independent columns give 1000 null experiments with one fixed design.
  null <- .pattern_mass_lm(matrix(rnorm(100 * 1000), 100), s)
  expect_lt(abs(mean(null$p_omnibus < 0.05) - 0.05), 0.025)
  expect_lt(abs(mean(null$p_omnibus) - 0.5), 0.035)
  power <- vapply(551:560, function(seed) {
    z <- confirmation_fixture(seed); c <- confirm_fixture(z)
    mean(c$omnibus$p_holm[c(1, 2, 3)] < 0.05)
  }, numeric(1))
  expect_gt(mean(power), 0.95)
})

test_that("single categorical component compares heads over independent blocks", {
  set.seed(556); n <- 80
  y <- factor(rep(c("a", "b"), 40))
  X <- cbind(as.numeric(y) + rnorm(n), rnorm(n), rnorm(n))
  fit <- .pattern_fit(X, y, rank = 1, control = pattern_control())
  test <- cbind(as.numeric(y) + rnorm(n), rnorm(n), rnorm(n))
  blocks <- rep(1:10, times = 1:10); blocks <- rep(blocks, length.out = n)
  c <- pattern_confirm(fit, test, y, blocks, confirmation_plan("sign_flip", 99),
    observation_ids = paste0("c", 1:n), discovery_ids = paste0("d", 1:n),
    feature_ids = paste0("v", 1:3), preprocessing_id = "BOLD", subject_id = "s1")
  result <- pattern_component_tests(fit, c, test, y,
    list(dataset = X, design = y, observation_ids = paste0("d", 1:n)))
  expect_equal(dim(result$losses$block_gain), c(10L, 1L))
  expect_equal(result$incremental$df, 9)
  expect_equal(result$incremental$gain, mean(vapply(split(seq_len(n), blocks), function(i)
    mean(result$losses$reduced[i, 1] - result$losses$full[i]), numeric(1))))
  expect_true(all(is.finite(result$association$p)))
  expect_true(all(is.finite(result$incremental$p)))
  expect_equal(c$prediction$metric, c("Accuracy", "logloss", "Brier"))
  expect_output(print(c), "3 features x 1 frozen")
  expect_error(pattern_component_tests(fit, c, test, y, list()), "calibration requires")
})

test_that("degenerate frozen decoding scores do not produce component significance", {
  d <- confirmation_fixture()
  d$fit$A[] <- d$fit$precision_A[] <- d$fit$G[] <- 0
  c <- confirm_fixture(d)
  out <- pattern_component_tests(d$fit, c, d$Xtest, d$Ytest)
  expect_true(all(is.na(out$association$p)))
  expect_error(pattern_component_tests(d$fit, c, d$Xtest, d$Ytest,
    list(dataset = d$X, design = d$Y, observation_ids = d$discovery)), "rank deficient")
})
