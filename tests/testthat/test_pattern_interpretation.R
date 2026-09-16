# Independent dense Gaussian references and adversarial locality cases.
phase4_quiet <- function(expr) { utils::capture.output(value <- suppressMessages(expr)); value }

phase4_fit <- function() {
  set.seed(440)
  X <- matrix(rnorm(120 * 6), 120, 6)
  Y <- cbind(X[, 1] + rnorm(120), X[, 2] + rnorm(120))
  .pattern_fit(X, Y, rank = 2, control = pattern_control(x_scale = "sd"))
}

rotate_fit_coordinates <- function(f, Q) {
  f$A <- f$A %*% Q; f$C <- f$C %*% Q
  f$precision_A <- f$precision_A %*% Q
  f$G <- crossprod(Q, f$G %*% Q)
  f$target_cov <- crossprod(Q, f$target_cov %*% Q)
  f
}

test_that("information agrees with dense log determinants including a suppressor", {
  f <- phase4_fit()
  f$A <- rbind(c(1, 0), c(0, 1), c(0, 0), c(0.5, 0.2), c(0, 0), c(0, 0))
  f$x_transform$sd <- c(2, 3, 4, 1, 1, 1)
  f$noise <- new_pattern_noise(rep(0.5, 6), matrix(c(2, 0, 2, 0, 0, 0), 6, 1))
  f$precision_A <- .noise_apply_precision(f$noise, f$A)
  f$G <- crossprod(f$A, f$precision_A)
  f$target_cov <- matrix(c(2, 0.4, 0.4, 1), 2)
  Psi <- diag(f$noise$D) + tcrossprod(f$noise$U)
  total <- Psi + f$A %*% f$target_cov %*% t(f$A)
  logdet <- function(M) as.numeric(determinant(M, logarithm = TRUE)$modulus)
  info <- (logdet(total) - logdet(Psi)) / 2
  reference <- vapply(1:6, function(v) info - (logdet(total[-v, -v]) - logdet(Psi[-v, -v])) / 2, numeric(1))
  expect_equal(model_importance(f, type = "conditional_info"), reference, tolerance = 1e-10)
  expect_equal(model_importance(f)[3], 0)
  expect_gt(reference[3], 0.1)
  expect_equal(model_importance(f), sqrt(diag(f$A %*% f$target_cov %*% t(f$A))) * f$x_transform$sd)
  Q <- matrix(c(cos(0.6), sin(0.6), -sin(0.6), cos(0.6)), 2)
  rotated <- rotate_fit_coordinates(f, Q)
  for (type in c("signal_sd", "conditional_info")) expect_equal(model_importance(f, type = type), model_importance(rotated, type = type), tolerance = 1e-12)
  # Singular target covariance still has a well-defined information map.
  f$target_cov <- diag(c(1, 0))
  expect_true(all(is.finite(model_importance(f, type = "conditional_info"))))
})

test_that("orthogonal views preserve fit predictions and invariant maps exactly", {
  f <- phase4_fit()
  set.seed(41); X <- matrix(rnorm(60), 10, 6)
  v <- rotate_patterns(f)
  expect_equal(v$L_b %*% v$H %*% t(v$L_t), tcrossprod(f$A, f$C), tolerance = 1e-12)
  expect_identical(predict(v, X, type = "decode"), predict(f, X, type = "decode"))
  expect_identical(predict(v, X, type = "scores"), predict(f, X, type = "scores"))
  for (type in c("signal_sd", "conditional_info")) expect_identical(model_importance(v, type = type), model_importance(f, type = type))
  expect_equal(model_patterns(v), model_patterns(f) %*% v$Q_b)
  Q <- matrix(c(0, 1, -1, 0), 2)
  rotated <- rotate_fit_coordinates(f, Q)
  expect_equal(predict(rotated, X, type = "decode"), predict(f, X, type = "decode"), tolerance = 1e-12)
  expect_equal(predict(rotated, X, type = "scores"), predict(f, X, type = "scores") %*% Q, tolerance = 1e-12)
  expect_error(rotate_patterns(f, spatial = diag(c(1, 2))), "orthogonal")
  expect_error(rotate_patterns(f, target = "oblimin"))
  different_targets <- f; different_targets$y_transform$mu <- f$y_transform$mu + 1
  expect_false(identical(rotate_patterns(different_targets)$basis_id, v$basis_id))
})

test_that("regional prediction uses marginal covariance and ignores outside columns", {
  f <- phase4_fit()
  f$noise <- new_pattern_noise(seq(0.5, 1, length.out = 6), matrix(c(1, 2, 1, 3, 2, 1), 6))
  f$precision_A <- .noise_apply_precision(f$noise, f$A); f$G <- crossprod(f$A, f$precision_A)
  region <- c(1L, 3L, 5L)
  local <- .pattern_restrict_fit(f, region)
  set.seed(42); X <- matrix(rnorm(60), 10, 6)
  Psi <- diag(f$noise$D) + tcrossprod(f$noise$U)
  PA <- solve(Psi[region, region], f$A[region, , drop = FALSE])
  G <- crossprod(f$A[region, ], PA)
  xc <- sweep(sweep(X[, region], 2, f$x_transform$mu[region], "-"), 2, f$x_transform$sd[region], "/")
  expect_equal(.pattern_sufficient(local, X), xc %*% PA, tolerance = 1e-10)
  That <- t(solve(solve(f$target_cov) + G, t(xc %*% PA)))
  expected <- .pattern_targets_invert(f$y_transform, That %*% t(f$C))
  expect_equal(predict(local, X, type = "decode"), expected, tolerance = 1e-10)
  Xbad <- X; Xbad[, -region] <- NA_real_
  expect_identical(predict(local, Xbad, type = "decode"), predict(local, X, type = "decode"))
  expect_gt(max(abs(PA - f$precision_A[region, ])), 1e-5)
  # A region containing only screened columns gives the training-mean predictor.
  f$p_input <- 7L
  empty <- .pattern_restrict_fit(f, 7L)
  pred <- predict(empty, cbind(X, 1), type = "decode")
  expect_equal(unname(pred), matrix(f$y_transform$mu, 10, 2, byrow = TRUE))
})

test_that("held-out Haufe matches a dense covariance oracle and checks overlap", {
  f <- phase4_fit()
  set.seed(43); X <- matrix(rnorm(600), 100, 6)
  f$training_observation_ids <- paste0("train:", 1:120)
  expect_error(pattern_haufe(f, X), "observation_ids")
  expect_error(pattern_haufe(f, X, paste0("train:", 1:100)), "overlap")
  h <- pattern_haufe(f, X, paste0("test:", 1:100))
  W <- (f$precision_A %*% solve(f$G)) / f$x_transform$sd
  S <- cov(X)
  oracle <- S %*% W %*% solve(crossprod(W, S %*% W))
  expect_equal(unname(h$empirical), unname(oracle), tolerance = 1e-10)
  Q <- matrix(c(0, 1, -1, 0), 2)
  hr <- pattern_haufe(rotate_fit_coordinates(f, Q), X, paste0("test:", 1:100))
  expect_equal(h$relative_discrepancy, hr$relative_discrepancy, tolerance = 1e-12)
})

test_that("subspace stability is rotation invariant and exposes unequal ranks", {
  f <- phase4_fit()
  Q <- matrix(c(0, 1, -1, 0), 2)
  g <- rotate_fit_coordinates(f, Q)
  result <- structure(list(fold_fits = list(f, g)), class = "pattern_global_result")
  s <- component_stability(result)
  expect_equal(s$overlap, 1, tolerance = 1e-12)
  expect_true(max(s$angles[[1]]) < 1e-7)
  g$A <- g$A[, 1, drop = FALSE]
  result$fold_fits[[2]] <- g
  expect_equal(component_stability(result)$overlap, 0.5, tolerance = 1e-12)
  g$A[] <- 0; result$fold_fits[[2]] <- g
  expect_true(is.na(component_stability(result)$overlap))
  result$fold_fits <- list(f, NULL, rotate_fit_coordinates(f, Q))
  skipped <- component_stability(result)
  expect_equal(skipped$fold1, 1L)
  expect_equal(skipped$fold2, 3L)
  expect_equal(skipped$overlap, 1, tolerance = 1e-12)
})

test_that("global regional ledgers, maps, diagnostics, and serialization agree", {
  sim <- sim_pattern_data(n = 60, dims = c(4, 4, 2), K = 3, seed = 444)
  res <- phase4_quiet(run_global(pattern_model(sim$dataset, sim$design, rank = 2), return_fits = TRUE, refit = TRUE))
  loc <- local_performance(res, list(all = 1:32, first = 1:8))
  expect_equal(loc$whole_brain[loc$region == "all"], loc$local_restricted[loc$region == "all"], tolerance = 1e-12)
  expect_equal(length(res$haufe_diagnostics), 3L)
  expect_equal(nrow(res$component_stability), 3L)
  expect_equal(res$fold_fits[[1]]$training_observation_ids,
               paste0("train:", setdiff(1:60, res$fold_ledger$observation[res$fold_ledger$fold == 1])))
  cmp <- local_performance(res, list(all = 1:32), independent_roi = list(all = res$fold_ledger))
  expect_equal(cmp$independent_roi, cmp$whole_brain)
  bad <- res$fold_ledger; bad$observation <- rev(bad$observation)
  expect_error(local_performance(res, list(all = 1:32), list(all = bad)), "must match")
  expect_error(local_performance(res, list(bad = c(1, 1))), "unique")
  expect_error(local_performance(res, list(bad = 0)), "valid")
  expect_error(local_performance(res, list(bad = 1e20)), "valid")
  invalid <- res$fold_ledger; invalid$prediction[] <- 2
  expect_error(local_performance(res, list(all = 1:32), list(all = invalid)), "probabilities")
  nofits <- res; nofits$fold_fits <- NULL
  expect_error(local_performance(nofits, list(all = 1:32)), "return_fits")
  map <- model_importance(res)
  expect_s4_class(map, "NeuroVol")
  expect_equal(as.numeric(map), as.numeric(model_importance(res$refit)))
  path <- tempfile("pattern-save-"); on.exit(unlink(path, recursive = TRUE))
  saved <- save_results(res, path, quiet = TRUE)
  restored <- readRDS(saved$result)
  expect_equal(restored$ledger, res$ledger)
  expect_identical(predict(restored$refit, sim$X, type = "prob"), predict(res$refit, sim$X, type = "prob"))
  expect_true(length(list.files(file.path(path, "maps"), pattern = "nii.gz")) == 2L)
})

test_that("continuous external assessment keeps training baselines under restriction", {
  sim <- sim_pattern_data(n = 45, dims = c(4, 4, 2), q = 2, r = 2,
                          external_test = TRUE, n_test = 21, seed = 445)
  res <- phase4_quiet(run_global(pattern_model(sim$dataset, sim$design, rank = 2), return_fits = TRUE))
  loc <- local_performance(res, list(all = 1:32, first = 1:8))
  expect_equal(loc$whole_brain[loc$region == "all"], loc$local_restricted[loc$region == "all"])
  ll <- attr(loc, "ledgers")$first
  expect_equal(ll$baseline, res$ledger$baseline)
  expect_equal(ll$partition, "external")
  expect_null(res$component_stability)
})

test_that("regional repeated CV uses one vote per observation", {
  sim <- sim_pattern_data(n = 60, dims = c(4, 4, 2), q = 2, r = 2, seed = 446)
  cv <- bootstrap_blocked_cross_validation(sim$design$block_var, nreps = 3)
  res <- phase4_quiet(run_global(pattern_model(sim$dataset, sim$design, rank = 2, crossval = cv), return_fits = TRUE))
  loc <- local_performance(res, list(first = 1:8, all = 1:32))
  expect_equal(loc$whole_brain[loc$region == "all"], loc$local_restricted[loc$region == "all"], tolerance = 1e-12)
  fl <- attr(loc, "fold_ledgers")$first
  pl <- attr(loc, "ledgers")$first
  expect_equal(pl$observation, sort(unique(fl$observation)))
  expected <- rowsum(fl$prediction, fl$observation) / as.numeric(table(fl$observation))
  expect_equal(unname(pl$prediction), unname(expected))
  expect_equal(pl$baseline, res$ledger$baseline)
})

test_that("screened columns remain missing in maps and rank-deficient scores are explicit", {
  set.seed(447)
  X <- cbind(matrix(rnorm(300), 60, 5), constant = 1)
  Y <- cbind(X[, 1], X[, 2])
  f <- .pattern_fit(X, Y, rank = 2)
  expect_true(is.na(model_importance(f)[6]))
  expect_true(all(is.na(model_patterns(f)[6, ])))
  f$A[, 2] <- 0
  f$precision_A <- .noise_apply_precision(f$noise, f$A)
  f$G <- crossprod(f$A, f$precision_A)
  h <- pattern_haufe(f, X)
  expect_equal(h$score_rank, 1)
  expect_true(all(is.finite(h$empirical)))
  expect_equal(unname(h$model[, 2]), rep(0, 5))
})

test_that("whole-domain Gaussian decoding has lower expected loss than restriction", {
  f <- phase4_fit()
  f$noise <- new_pattern_noise(rep(0.5, 6), matrix(c(2, 0, 2, 0, 0, 0), 6))
  f$precision_A <- .noise_apply_precision(f$noise, f$A)
  f$G <- crossprod(f$A, f$precision_A)
  region <- .pattern_restrict_fit(f, 1:2)
  risk <- function(G) sum(diag(solve(solve(f$target_cov) + G)))
  expect_lt(risk(f$G), risk(region$G))
  set.seed(448)
  T <- matrix(rnorm(20000 * 2), 20000, 2)
  E <- matrix(rnorm(20000 * 6), 20000, 6) * sqrt(0.5) + tcrossprod(rnorm(20000), f$noise$U[, 1])
  Xc <- T %*% t(f$A) + E
  full <- Xc %*% f$precision_A %*% solve(diag(2) + f$G)
  local <- Xc[, 1:2] %*% region$precision_A %*% solve(diag(2) + region$G)
  expect_lt(mean((full - T)^2), mean((local - T)^2))
})


test_that("feature units scale forward maps but not conditional information", {
  set.seed(449)
  X <- matrix(rnorm(100 * 6), 100, 6)
  Y <- cbind(X[, 1] + X[, 2], X[, 3])
  units <- c(0.01, 2, 100, 4, 0.1, 3)
  ctl <- pattern_control(x_scale = "sd", noise = list(type = "identity"))
  a <- .pattern_fit(X, Y, rank = 2, control = ctl)
  b <- .pattern_fit(sweep(X, 2, units, "*"), Y, rank = 2, control = ctl)
  expect_equal(model_importance(b), model_importance(a) * units, tolerance = 1e-10)
  expect_equal(model_importance(b, type = "conditional_info"),
               model_importance(a, type = "conditional_info"), tolerance = 1e-10)
})


test_that("Haufe score rank uses the same spectral cutoff as calibrated prediction", {
  f <- phase4_fit()
  Q <- matrix(c(1, 1, -1, 1), 2) / sqrt(2)
  f$noise <- new_pattern_noise(rep(1, 6), type = "identity")
  f$A <- rbind(diag(sqrt(c(1, 7e-11))), matrix(0, 4, 2)) %*% t(Q)
  f$precision_A <- f$A; f$G <- crossprod(f$A)
  set.seed(450); X <- matrix(rnorm(600), 100, 6)
  expect_equal(pattern_haufe(f, X)$score_rank, 1L)
})
