library(testthat)

quiet_run <- function(expr) {
  utils::capture.output(res <- suppressMessages(expr))
  res
}

test_that("the group-lasso prox and lambda_max behave as advertised", {
  set.seed(1)
  B <- matrix(rnorm(20 * 3), 20, 3)
  nrm <- sqrt(rowSums(B^2))
  thr <- stats::median(nrm)
  P <- rMVPA:::.prox_group_lasso(B, thr)
  # rows shorter than the threshold are zeroed exactly; the rest shrink toward 0
  expect_true(all(rowSums(P[nrm <= thr, , drop = FALSE]^2) == 0))
  kept <- which(nrm > thr)
  expect_equal(sqrt(rowSums(P[kept, , drop = FALSE]^2)), nrm[kept] - thr, tolerance = 1e-12)
  # direction is preserved
  expect_equal(P[kept[1], ] / sqrt(sum(P[kept[1], ]^2)),
               B[kept[1], ] / nrm[kept[1]], tolerance = 1e-12)
  expect_equal(rMVPA:::.prox_group_lasso(B, 0), B)

  # lambda_max is by definition the smallest penalty that zeroes every feature
  # at A = 0, which is a statement about one prox step, so test it there.
  X <- matrix(rnorm(40 * 15), 40, 15)
  y <- factor(rep(c("a", "b"), 20))
  fit <- rMVPA:::.pattern_fit(X, y, rank = 1)
  Xc <- scale(X, scale = FALSE)
  Tm <- rMVPA:::.pattern_targets_apply(fit$y_transform, y) %*% fit$C
  const <- rMVPA:::.noise_quadform(fit$noise, Xc)
  pen <- rMVPA:::.pattern_resolve_penalty(list(sparse = 0.5), Xc, Tm, fit$noise, NULL)
  expect_true(is.finite(pen$lambda_max))
  expect_equal(pen$lambda_s, 0.5 * pen$lambda_max)
  expect_error(rMVPA:::.pattern_resolve_penalty(list(sparse = 1), Xc, Tm, fit$noise, NULL),
               "\\[0, 1\\)")
  A0 <- matrix(0, ncol(Xc), 1)
  step_at <- function(lambda_s) {
    prob <- rMVPA:::.pattern_astep_problem(Xc, Tm, fit$noise, NULL, lambda_s, 0, const)
    st <- 1 / prob$lipschitz
    rMVPA:::.prox_group_lasso(A0 - st * rMVPA:::.pattern_astep_gradient(A0, prob), st * lambda_s)
  }
  expect_true(all(step_at(pen$lambda_max * (1 + 1e-9)) == 0))       # at lambda_max: all zero
  expect_gt(sum(rowSums(step_at(pen$lambda_max * 0.99)^2) > 0), 0L) # just below: something survives

  # and the fitted support shrinks as the penalty grows
  sizes <- vapply(c(0.05, 0.2, 0.5), function(a) {
    sum(rowSums(rMVPA:::.pattern_fit(X, y, rank = 1, penalty = list(sparse = a))$A^2) > 0)
  }, numeric(1))
  expect_true(all(diff(sizes) <= 0))
  expect_gt(sizes[1], sizes[3])
})

test_that("the penalized objective decreases monotonically, inside and across steps", {
  sim <- sim_pattern_data(n = 90, dims = c(8, 8, 4), K = 3, snr = 0.6, seed = 60)
  g <- spatial_graph(sim$dataset)
  for (pen in list(list(sparse = 0.3),
                   list(signed_smooth = 1),
                   list(sparse = 0.3, signed_smooth = 0.5))) {
    fit <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 2, graph = g, penalty = pen)
    tr <- fit$diagnostics$objective
    info <- paste(names(pen), unlist(pen), sep = "=", collapse = ",")
    expect_true(all(diff(tr) <= 1e-9 * abs(tr[1])), info = info)   # outer alternation
    expect_true(fit$diagnostics$converged, info = info)
    expect_gt(length(fit$diagnostics$inner), 0L)
  }

  # the inner solver's own trace is non-increasing (monotone FISTA)
  Xc <- scale(sim$X, scale = FALSE)
  fit <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 2, graph = g,
                              penalty = list(sparse = 0.2, signed_smooth = 0.5))
  Yw <- rMVPA:::.pattern_targets_apply(fit$y_transform, sim$targets)
  prob <- rMVPA:::.pattern_astep_problem(Xc, Yw %*% fit$C, fit$noise, fit$penalty$L,
                                         fit$penalty$lambda_s, fit$penalty$lambda_l,
                                         rMVPA:::.noise_quadform(fit$noise, Xc))
  sol <- rMVPA:::.pattern_astep_fista(matrix(0, ncol(Xc), 2), prob, max_iter = 100)
  expect_true(all(diff(sol$trace) <= 1e-12 * abs(sol$trace[1])))
  expect_true(sol$converged)

  # the cheap A-step objective agrees with the explicit residual form
  A <- sol$A
  ref <- rMVPA:::.pattern_objective(Xc, Yw, A, fit$C, fit$noise, fit$penalty)
  expect_equal(rMVPA:::.pattern_astep_objective(A, prob), ref, tolerance = 1e-8)
})

test_that("the A-step gradient matches a finite-difference reference", {
  set.seed(2)
  n <- 40; p <- 25; r <- 2
  X <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  Tm <- matrix(rnorm(n * r), n, r)
  noise <- rMVPA:::new_pattern_noise(runif(p, 0.5, 2), matrix(rnorm(p * 2), p, 2) * 0.5)
  A0 <- as.logical(diag(p)[, 1])
  g <- spatial_graph(list(A = as.matrix(Matrix::bandSparse(p, k = 1, diagonals = list(rep(1, p - 1))))))
  prob <- rMVPA:::.pattern_astep_problem(X, Tm, noise, g$L / 4, 0, 0.7,
                                         rMVPA:::.noise_quadform(noise, X))
  A <- matrix(rnorm(p * r), p, r)
  G <- rMVPA:::.pattern_astep_gradient(A, prob)
  # only the smooth part has a gradient, so compare with lambda_s = 0
  smooth_obj <- function(M) rMVPA:::.pattern_astep_objective(M, prob)
  eps <- 1e-6
  fd <- matrix(0, p, r)
  for (i in c(1, 7, p)) for (j in seq_len(r)) {
    Ap <- A; Am <- A
    Ap[i, j] <- Ap[i, j] + eps; Am[i, j] <- Am[i, j] - eps
    fd[i, j] <- (smooth_obj(Ap) - smooth_obj(Am)) / (2 * eps)
  }
  for (i in c(1, 7, p)) for (j in seq_len(r)) {
    expect_equal(G[i, j], fd[i, j], tolerance = 1e-5)
  }
})

test_that("sparsity recovers the informative support on known-truth scenarios", {
  # precision = share of selected features that are truly informative
  # recall    = share of informative features that were selected
  floors <- list(basic = c(0.9, 0.6), redundant = c(0.9, 0.6))
  for (sc in names(floors)) {
    sim <- sim_pattern_data(n = 120, dims = c(8, 8, 4), K = 3, snr = 0.8,
                            scenario = sc, seed = 61)
    fit <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = sim$r,
                                penalty = list(sparse = 0.1))
    sel <- which(rowSums(fit$A^2) > 0)
    precision <- mean(sel %in% sim$informative)
    recall <- mean(sim$informative %in% sel)
    expect_gt(length(sel), 0L)
    expect_lt(length(sel), nrow(fit$A))          # something was actually excluded
    expect_gt(precision, floors[[sc]][1], label = paste(sc, "precision"))
    expect_gt(recall, floors[[sc]][2], label = paste(sc, "recall"))
  }
  # the redundant scenario expresses one dimension in two separate regions;
  # both must be found, not just the larger one
  sim <- sim_pattern_data(n = 120, dims = c(8, 8, 4), K = 3, snr = 0.8,
                          scenario = "redundant", seed = 61)
  fit <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 1, penalty = list(sparse = 0.1))
  sel <- which(rowSums(fit$A^2) > 0)
  expect_gt(mean(intersect(sim$boxes[[1]], sim$informative) %in% sel), 0.5)
  expect_gt(mean(intersect(sim$boxes[[4]], sim$informative) %in% sel), 0.5)
})

test_that("the suppressor region is excluded from the patterns", {
  # x1 = s + n, x2 = n: the second region cancels noise but expresses no task
  # signal, so a forward-pattern penalty should leave it out.
  sim <- sim_pattern_data(n = 120, dims = c(8, 8, 4), K = 3, snr = 0.8,
                          scenario = "suppressor", seed = 62)
  fit <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 1, penalty = list(sparse = 0.25))
  sel <- which(rowSums(fit$A^2) > 0)
  suppressor_only <- setdiff(sim$boxes[[2]], sim$informative)
  expect_gt(length(suppressor_only), 0L)
  # far fewer suppressor voxels survive than informative ones
  expect_lt(mean(suppressor_only %in% sel), mean(sim$informative %in% sel))
})

test_that("signed smoothing helps a smooth pattern and erases a sign-flipping one", {
  # The plan's key check: 'sparse and smooth' must not silently become
  # 'incapable of representing a spatially heterogeneous code'.
  # signed_smooth is measured against the data-fit curvature, so these are
  # dimensionless multiples of "as influential as the data".
  study <- function(scenario, rhos = c(0, 4, 16)) {
    sim <- sim_pattern_data(n = 120, dims = c(8, 8, 4), K = 3, snr = 0.8,
                            scenario = scenario, seed = 63)
    g <- spatial_graph(sim$dataset)
    tr <- seq_len(80); te <- 81:120
    vapply(rhos, function(rho) {
      fit <- rMVPA:::.pattern_fit(sim$X[tr, ], sim$targets[tr], rank = sim$r, graph = g,
                                  penalty = list(sparse = 0.15, signed_smooth = rho))
      c(loss = rMVPA:::.pattern_loss(fit, sim$X[te, ], sim$targets[te]),
        pattern_cor = abs(stats::cor(as.numeric(fit$A[, 1]), as.numeric(sim$A[, 1]))))
    }, numeric(2))
  }
  smooth <- study("smooth")     # a coherent blob: smoothing is the right prior
  flip <- study("signflip")     # a contiguous region of alternating signs

  # smoothing recovers a coherent blob better
  expect_gt(smooth["pattern_cor", 2], smooth["pattern_cor", 1] + 0.05)
  # and degrades a fine-grained sign-flipping code monotonically
  expect_lt(flip["pattern_cor", 3], flip["pattern_cor", 1] - 0.1)
  expect_true(all(diff(flip["pattern_cor", ]) < 0))
  expect_gt(flip["pattern_cor", 1], 0.9)
  expect_lte(flip["loss", 1], flip["loss", 3])
  # Honest caveat worth pinning: on the smooth scenario smoothing improves the
  # recovered pattern while making held-out decoding worse. Better anatomy and
  # better prediction are not the same objective, so both are reported.
  expect_gt(smooth["loss", 2], smooth["loss", 1])
  # the two scenarios respond to the same smoothing in opposite directions:
  # this is why signed smoothing is tunable and defaults to off
  expect_gt(smooth["pattern_cor", 2] - smooth["pattern_cor", 1],
            flip["pattern_cor", 2] - flip["pattern_cor", 1])
})

test_that("the penalized fit is reproducible", {
  sim <- sim_pattern_data(n = 60, dims = c(5, 5, 3), K = 3, snr = 0.8, seed = 68)
  g <- spatial_graph(sim$dataset)
  pen <- list(sparse = 0.2, signed_smooth = 0.5)
  a <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 2, graph = g, penalty = pen)
  b <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 2, graph = g, penalty = pen)
  expect_identical(a$A, b$A)
  expect_identical(a$diagnostics$objective, b$diagnostics$objective)
  # and it does not disturb the global RNG stream
  set.seed(1); before <- stats::runif(1)
  set.seed(1); invisible(rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 2, graph = g, penalty = pen))
  expect_equal(stats::runif(1), before)
})

test_that("permuting features with the graph leaves the penalized fit unchanged", {
  sim <- sim_pattern_data(n = 90, dims = c(6, 6, 3), K = 3, snr = 0.8, seed = 64)
  g <- spatial_graph(sim$dataset)
  pen <- list(sparse = 0.2, signed_smooth = 0.5)
  fit <- rMVPA:::.pattern_fit(sim$X, sim$targets, rank = 2, graph = g, penalty = pen)

  set.seed(9)
  perm <- sample(ncol(sim$X))
  gp <- restrict_graph(g, perm)
  fitp <- rMVPA:::.pattern_fit(sim$X[, perm], sim$targets, rank = 2, graph = gp, penalty = pen)

  # The problem is permutation equivariant, and the fit is deterministic, but
  # the step-size power iteration starts from a fixed vector that is not itself
  # permutation equivariant, so the two runs stop at very slightly different
  # points on the same solution path.
  expect_equal(fitp$A, fit$A[perm, , drop = FALSE], tolerance = 1e-4)
  expect_equal(fitp$C, fit$C, tolerance = 1e-3)
  expect_equal(fitp$penalty$lambda_s, fit$penalty$lambda_s, tolerance = 1e-8)
  expect_equal(tail(fitp$diagnostics$objective, 1), tail(fit$diagnostics$objective, 1),
               tolerance = 1e-3)
  # the selected support is identical
  expect_equal(which(rowSums(fitp$A^2) > 0), which(rowSums(fit$A[perm, , drop = FALSE]^2) > 0))
})

test_that("the graph is restricted when column screening drops features", {
  sim <- sim_pattern_data(n = 60, dims = c(5, 5, 3), K = 2, snr = 1, seed = 65)
  X <- sim$X
  X[, c(4, 20)] <- 5                                  # constant columns
  g <- spatial_graph(sim$dataset)
  fit <- rMVPA:::.pattern_fit(X, sim$targets, rank = 1, graph = g,
                              penalty = list(sparse = 0.2, signed_smooth = 0.5))
  expect_equal(fit$feature_index, setdiff(seq_len(ncol(X)), c(4, 20)))
  expect_equal(nrow(fit$A), ncol(X) - 2L)
  expect_equal(fit$graph$n_features, ncol(X) - 2L)
  expect_equal(fit$graph$feature_ids, g$feature_ids[fit$feature_index])
  # a mismatched graph is refused rather than silently misaligned
  expect_error(rMVPA:::.pattern_fit(X[, 1:10], sim$targets, rank = 1,
                                    graph = restrict_graph(g, 1:12),
                                    penalty = list(signed_smooth = 1)),
               "12 vertices")
})

test_that("smoothing without a graph is refused", {
  set.seed(3)
  X <- matrix(rnorm(40 * 12), 40, 12)
  y <- factor(rep(c("a", "b"), 20))
  expect_error(rMVPA:::.pattern_fit(X, y, rank = 1, penalty = list(signed_smooth = 1)),
               "requires a spatial graph")
  # sparsity alone needs no anatomy
  expect_s3_class(rMVPA:::.pattern_fit(X, y, rank = 1, penalty = list(sparse = 0.3)),
                  "pattern_fit")
})

test_that("penalty selection is nested and no assessment row reaches the fit", {
  sim <- sim_pattern_data(n = 90, dims = c(6, 6, 3), K = 3, snr = 0.8, seed = 66)
  spec <- pattern_model(sim$dataset, sim$design, max_rank = 2,
                        penalty = list(sparse = "auto"), return_fits = TRUE)
  res <- quiet_run(run_global(spec, return_fits = TRUE))
  expect_true(all(res$alphas > 0))
  expect_true(all(res$n_nonzero > 0 & res$n_nonzero <= ncol(get_feature_matrix(sim$dataset))))
  expect_true("n_selected" %in% names(performance(res)))
  expect_equal(names(output_schema(spec)),
               c("Accuracy", "AUC", "logloss", "rank_mean", "n_selected"))

  # poison every fold-1 test row: the fold-1 fit must be untouched
  poisoned <- sim
  te <- which(sim$design$block_var == 1)
  Xp <- sim$X; Xp[te, ] <- Xp[te, ] * 1000 + 777
  arr <- array(t(Xp), c(sim$dims, nrow(Xp)))
  poisoned$dataset <- mvpa_dataset(
    train_data = neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(sim$dims, nrow(Xp)), c(1, 1, 1))),
    mask = sim$dataset$mask
  )
  spec_p <- pattern_model(poisoned$dataset, sim$design, max_rank = 2,
                          penalty = list(sparse = "auto"), return_fits = TRUE)
  res_p <- quiet_run(run_global(spec_p, return_fits = TRUE))
  f1 <- res$fold_fits[[1]]; f1p <- res_p$fold_fits[[1]]
  # identical training rows must give a bit-identical fit, penalty scale included
  expect_identical(f1p$A, f1$A)
  expect_identical(f1p$penalty$lambda_s, f1$penalty$lambda_s)
  expect_identical(f1p$penalty$lambda_max, f1$penalty$lambda_max)
  expect_identical(f1p$x_transform$mu, f1$x_transform$mu)
  expect_identical(f1p$noise$D, f1$noise$D)
})

test_that("a penalized model runs end to end in regional and searchlight modes", {
  sim <- sim_pattern_data(n = 60, dims = c(6, 6, 3), K = 3, snr = 1, seed = 67)
  region <- neuroim2::NeuroVol(array(rep(1:3, length.out = 108), sim$dims),
                               neuroim2::space(sim$dataset$mask))
  spec <- pattern_model(sim$dataset, sim$design, rank = 1,
                        penalty = list(sparse = 0.2, signed_smooth = 0.5))
  expect_s3_class(spec$graph, "spatial_graph")
  reg <- quiet_run(run_regional(spec, region))
  expect_equal(nrow(reg$performance_table), 3L)
  expect_true("n_selected" %in% names(reg$performance_table))
  expect_true(all(reg$performance_table$n_selected > 0))
  # each ROI selects from its own features only
  expect_true(all(reg$performance_table$n_selected <= 36))

  sl <- quiet_run(run_searchlight(spec, radius = 2, method = "standard"))
  expect_true("n_selected" %in% names(sl$results))
})

test_that("the graph stays sparse at whole-brain feature counts", {
  # Symmetrizing with base::pmax would coerce to a dense p x p matrix and
  # exhaust memory here; the sparse identity max(a,b) = (a+b+|a-b|)/2 must be used.
  ds <- gen_sample_dataset(c(40, 40, 20), 4)
  t0 <- proc.time()[["elapsed"]]
  g <- spatial_graph(ds$dataset)
  elapsed <- proc.time()[["elapsed"]] - t0
  expect_equal(g$n_features, 32000L)
  expect_s4_class(g$A, "dgCMatrix")
  expect_lt(as.numeric(object.size(g$A)) / 1024^2, 50)   # far below a dense 32000^2
  expect_lt(elapsed, 30)
  expect_equal(g$degree, as.numeric(Matrix::rowSums(g$A)))
  expect_true(is.finite(g$L_norm) && g$L_norm > 0)
})

test_that("signed_smooth is scaled to the data and does not drift with sample size", {
  sim <- sim_pattern_data(n = 120, dims = c(6, 6, 3), K = 3, snr = 0.8, seed = 70)
  g <- spatial_graph(sim$dataset)
  for (nn in c(60, 120)) {
    fit <- rMVPA:::.pattern_fit(sim$X[seq_len(nn), ], sim$targets[seq_len(nn)], rank = 1,
                                graph = g, penalty = list(signed_smooth = 1))
    Xc <- scale(sim$X[seq_len(nn), ], scale = FALSE)
    Yw <- rMVPA:::.pattern_targets_apply(fit$y_transform, sim$targets[seq_len(nn)])
    curvature <- rMVPA:::.pattern_grad_norm(fit$noise, crossprod(Yw %*% fit$C), NULL, 0,
                                            nn, ncol(Xc), 1L, safety = 1)
    # rho = 1 means "as influential as the data fit", at any sample size
    expect_equal(fit$penalty$lambda_l / curvature, 1, tolerance = 1e-6,
                 label = paste("n =", nn))
  }
})

test_that("the A-step step size tracks the true Lipschitz constant", {
  set.seed(71)
  n <- 60; p <- 80; r <- 2
  X <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  Tm <- matrix(rnorm(n * r), n, r)
  noise <- rMVPA:::new_pattern_noise(runif(p, 0.5, 2), matrix(rnorm(p * 2), p, 2) * 0.5)
  prob <- rMVPA:::.pattern_astep_problem(X, Tm, noise, NULL, 0.1, 0,
                                         rMVPA:::.noise_quadform(noise, X))
  # build the gradient operator explicitly and take its true spectral norm
  M <- matrix(0, p * r, p * r)
  for (j in seq_len(p * r)) {
    E <- matrix(0, p, r); E[j] <- 1
    M[, j] <- as.numeric(rMVPA:::.noise_apply_precision(noise, E %*% prob$TtT) / n)
  }
  true_L <- max(svd(M)$d)
  # an underestimate would break the descent guarantee, so it must be >= 1
  expect_gte(prob$lipschitz / true_L, 1)
  expect_lt(prob$lipschitz / true_L, 1.1)
})

test_that("the A-step reaches the same optimum as a slow reference solver", {
  set.seed(72)
  n <- 80; p <- 120; r <- 2
  X <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  X[, 1:20] <- X[, 1:20] + matrix(rnorm(n), n, 1) %*% matrix(1, 1, 20)
  Tm <- matrix(rnorm(n * r), n, r)
  noise <- rMVPA:::new_pattern_noise(runif(p, 0.5, 2), NULL)
  prob <- rMVPA:::.pattern_astep_problem(X, Tm, noise, NULL, 0.05, 0,
                                         rMVPA:::.noise_quadform(noise, X))
  # plain proximal gradient with a definitely-safe step, run to death
  step <- 1 / (2 * prob$lipschitz)
  Aref <- matrix(0, p, r)
  for (i in seq_len(20000)) {
    Aref <- rMVPA:::.prox_group_lasso(
      Aref - step * rMVPA:::.pattern_astep_gradient(Aref, prob), step * prob$lambda_s)
  }
  sol <- rMVPA:::.pattern_astep_fista(matrix(0, p, r), prob, max_iter = 500,
                                      tol = 1e-9, tol_iterate = 1e-6)
  expect_equal(sol$objective, rMVPA:::.pattern_astep_objective(Aref, prob), tolerance = 1e-9)
  expect_lt(sqrt(sum((sol$A - Aref)^2)) / sqrt(sum(Aref^2)), 1e-4)
  expect_equal(which(rowSums(sol$A^2) > 0), which(rowSums(Aref^2) > 0))
  expect_lt(sol$iterations, 200L)
})

test_that("a step that is too long is recovered from, not reported as convergence", {
  set.seed(73)
  n <- 60; p <- 80; r <- 2
  X <- scale(matrix(rnorm(n * p), n, p), scale = FALSE)
  Tm <- matrix(rnorm(n * r), n, r)
  noise <- rMVPA:::new_pattern_noise(runif(p, 0.5, 2), NULL)
  prob <- rMVPA:::.pattern_astep_problem(X, Tm, noise, NULL, 0.05, 0,
                                         rMVPA:::.noise_quadform(noise, X))
  good <- rMVPA:::.pattern_astep_fista(matrix(0, p, r), prob, max_iter = 500,
                                       tol = 1e-9, tol_iterate = 1e-6)
  bad <- prob; bad$lipschitz <- prob$lipschitz / 4       # step four times too long
  sol <- rMVPA:::.pattern_astep_fista(matrix(0, p, r), bad, max_iter = 500,
                                      tol = 1e-9, tol_iterate = 1e-6)
  # it must not stop immediately claiming success at the starting point
  expect_gt(sol$iterations, 1L)
  expect_equal(sol$objective, good$objective, tolerance = 1e-6)
  expect_gt(sum(rowSums(sol$A^2) > 0), 0L)
  expect_true(all(diff(sol$trace) <= 1e-12 * abs(sol$trace[1])))
})

test_that("multibasis ROIs map onto the right basis channel of the graph", {
  ds <- gen_sample_dataset(c(5, 5, 3), 40, nlevels = 2, blocks = 2)
  mb <- mvpa_multibasis_dataset(
    train_data = list(ds$dataset$train_data, ds$dataset$train_data),
    mask = ds$dataset$mask
  )
  g <- spatial_graph(mb)
  n_vox <- sum(ds$dataset$mask > 0)
  expect_true(anyDuplicated(g$feature_ids) > 0)          # ids repeat across channels

  # an ROI covering both channels of the same voxels must map to both channels
  ids <- rep(c(1L, 2L, 3L), times = 2)
  pos <- rMVPA:::.pattern_graph_positions(g, ids)
  expect_equal(pos, c(1L, 2L, 3L, n_vox + 1L, n_vox + 2L, n_vox + 3L))
  expect_false(anyDuplicated(pos) > 0)
  expect_equal(g$basis[pos], rep(1:2, each = 3))

  # and a regional run with smoothing completes instead of erroring every ROI
  spec <- pattern_model(mb, ds$design, rank = 1,
                        penalty = list(sparse = 0.2, signed_smooth = 1))
  region <- neuroim2::NeuroVol(array(rep(1:2, length.out = 75), c(5, 5, 3)),
                               neuroim2::space(ds$dataset$mask))
  reg <- quiet_run(run_regional(spec, region))
  expect_equal(nrow(reg$performance_table), 2L)
  expect_true(all(is.finite(reg$performance_table$Accuracy)))
})

test_that("tuning is deterministic even when the design has no blocks", {
  sim <- sim_pattern_data(n = 60, dims = c(5, 5, 3), K = 3, snr = 0.8, seed = 74)
  des <- mvpa_design(data.frame(y = sim$targets), y_train = ~ y)   # no block_var
  spec <- suppressWarnings(pattern_model(sim$dataset, des, max_rank = 2,
                                         penalty = list(sparse = "auto")))
  X <- get_feature_matrix(sim$dataset)
  a <- rMVPA:::.pattern_select_config(X, sim$targets, NULL, spec$control,
                                      penalty = spec$penalty)
  b <- rMVPA:::.pattern_select_config(X, sim$targets, NULL, spec$control,
                                      penalty = spec$penalty)
  expect_identical(a$rank, b$rank)
  expect_identical(a$penalty, b$penalty)
  expect_equal(a$losses, b$losses)
  # and it leaves the global RNG stream alone
  set.seed(1); before <- stats::runif(1)
  set.seed(1)
  invisible(rMVPA:::.pattern_select_config(X, sim$targets, NULL, spec$control,
                                           penalty = spec$penalty))
  expect_equal(stats::runif(1), before)
})
