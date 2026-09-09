# Synthetic subject estimates provide a direct oracle independent of the
# confirmation regression; full covariance is intentionally non-diagonal.
group_fixture <- function(s = 12, r = 2, p = 3) {
  set.seed(660)
  B <- diag(r)
  basis <- structure(list(matrix = B, target_ids = paste0("y", 1:r), type = "continuous", basis_id = "fixed"), class = "pattern_basis")
  lapply(seq_len(s), function(j) {
    A <- matrix(rnorm(p * r, sd = 0.4), p, r) + 0.6
    V <- diag(seq_len(r) / 20, nrow = r) + matrix(0.02, r, r)
    structure(list(estimate = A, covariance = array(rep(V, p), c(r, r, p)),
      basis = basis, feature_ids = paste0("v", 1:p), nuisance = matrix(1, 10, 1, dimnames = list(NULL, "intercept")),
      provenance = list(subject_id = paste0("s", j), observation_ids = paste0("s", j, ":c:", 1:10),
        discovery_ids = paste0("s", j, ":d:", 1:10), preprocessing_id = "BOLD")), class = "pattern_confirmation")
  })
}

test_that("random group estimates match equal-subject and Hotelling oracles", {
  ss <- group_fixture(); g <- pattern_group(ss, ss[[1]]$basis)
  Y <- do.call(rbind, lapply(ss, function(x) x$estimate[1, ]))
  V <- cov(Y) / 12; mu <- colMeans(Y)
  T2 <- as.numeric(t(mu) %*% solve(V) %*% mu)
  expect_equal(unname(g$estimate[1, ]), mu)
  expect_equal(g$covariance[, , 1], V)
  expect_equal(g$omnibus$statistic[1], T2 * 10/(2*11))
  expect_equal(g$omnibus$p[1], pf(T2 * 10/(2*11), 2, 10, lower.tail = FALSE))
  expect_equal(unname(g$p[1, 1]), t.test(Y[, 1])$p.value)
  expect_equal(g$heterogeneity$trace[1], sum(pmax(eigen(cov(Y) - ss[[1]]$covariance[, , 1], symmetric = TRUE)$values, 0)))
  expect_equal(g$subject_expression[1, 1], sum(Y[1, ] * colMeans(Y[-1, ])) / sqrt(sum(colMeans(Y[-1, ])^2)), ignore_attr = TRUE)
})

test_that("fixed group GLS preserves cross-component covariance", {
  ss <- group_fixture()
  for (j in seq_along(ss)) ss[[j]]$covariance <- ss[[j]]$covariance * j
  g <- pattern_group(ss, ss[[1]]$basis, effects = "fixed")
  Y <- unlist(lapply(ss, function(x) x$estimate[1, ]))
  V <- as.matrix(Matrix::bdiag(lapply(ss, function(x) x$covariance[, , 1])))
  D <- kronecker(matrix(1, 12, 1), diag(2))
  oracle_cov <- solve(t(D) %*% solve(V, D))
  oracle <- oracle_cov %*% t(D) %*% solve(V, Y)
  expect_equal(unname(g$estimate[1, ]), as.vector(oracle), tolerance = 1e-10)
  expect_equal(g$covariance[, , 1], oracle_cov, tolerance = 1e-10)
  expect_equal(g$heterogeneity$Q[1], as.numeric(crossprod(Y - D %*% oracle, solve(V, Y - D %*% oracle))), tolerance = 1e-10)
})

test_that("basis transport and orthogonal invariants survive heterogeneous coordinates", {
  ss <- group_fixture(); ref <- ss[[1]]$basis
  g <- pattern_group(ss, ref)
  moved <- ss
  for (j in seq_along(ss)) {
    H <- matrix(c(1 + j/10, 0.1, 0.3, 2), 2)
    moved[[j]]$basis$matrix <- ref$matrix %*% H
    moved[[j]]$estimate <- ss[[j]]$estimate %*% solve(t(H))
    for (v in 1:3) moved[[j]]$covariance[, , v] <- solve(H) %*% ss[[j]]$covariance[, , v] %*% solve(t(H))
  }
  aligned <- pattern_group(moved, ref)
  expect_equal(aligned$estimate, g$estimate, tolerance = 1e-10)
  expect_equal(aligned$covariance, g$covariance, tolerance = 1e-10)
  Q <- matrix(c(cos(.7), sin(.7), -sin(.7), cos(.7)), 2)
  rotated <- ref; rotated$matrix <- ref$matrix %*% Q
  for (mode in c("fixed", "random")) {
    original <- pattern_group(ss, ref, effects = mode)
    gr <- pattern_group(ss, rotated, effects = mode)
    expect_equal(unname(gr$estimate), unname(original$estimate %*% Q), tolerance = 1e-10)
    expect_equal(gr$omnibus, original$omnibus, tolerance = 1e-10)
    expect_equal(gr$effect_norm, original$effect_norm, tolerance = 1e-10)
    expect_equal(gr$heterogeneity, original$heterogeneity, tolerance = 1e-10)
    expect_equal(gr$subject_expression, original$subject_expression, tolerance = 1e-10)
  }
})

test_that("spatial correspondence is identity-safe and excludes unsupported interpolation", {
  ss <- group_fixture(); ref <- ss[[1]]$basis
  original <- pattern_group(ss, ref)
  moved <- ss
  for (j in seq_along(ss)) {
    moved[[j]]$estimate <- ss[[j]]$estimate[3:1, ]
    moved[[j]]$covariance <- ss[[j]]$covariance[, , 3:1]
    moved[[j]]$feature_ids <- ss[[j]]$feature_ids[3:1]
  }
  g <- pattern_group(moved, ref)
  expect_equal(g$estimate[3:1, ], original$estimate)
  mappings <- lapply(ss, function(x) setNames(x$feature_ids, paste0("shared", 1:3)))
  expect_equal(unname(pattern_group(ss, ref, mappings)$estimate), unname(original$estimate))
  mappings[[1]][2] <- mappings[[1]][1]
  expect_error(pattern_group(ss, ref, mappings), "one-to-one")
  bad <- ss; bad[[1]]$provenance$subject_id <- "s2"
  expect_error(pattern_group(bad, ref), "unique")
  bad <- ss; bad[[1]]$provenance$discovery_ids[1] <- bad[[2]]$provenance$observation_ids[1]
  expect_error(pattern_group(bad, ref), "overlap")
  bad <- ss; bad[[1]]$provenance$preprocessing_id <- "different-units"
  expect_error(pattern_group(bad, ref), "preprocessing")
  bad <- ss; bad[[1]]$basis$matrix <- matrix(1, 2, 2)
  expect_error(pattern_group(bad, ref), "full-rank")
  expect_error(pattern_group(ss[1:3], ref), "more subjects")
})

test_that("one-dimensional and singular group inference have honest limits", {
  ss <- group_fixture(r = 1, p = 1)
  g <- pattern_group(ss, ss[[1]]$basis)
  Y <- vapply(ss, function(x) x$estimate[1], numeric(1))
  expect_equal(g$omnibus$p, t.test(Y)$p.value)
  for (j in seq_along(ss)) ss[[j]]$estimate[] <- 1
  g <- pattern_group(ss, ss[[1]]$basis)
  expect_true(is.na(g$omnibus$p))
  expect_true(is.na(g$p[1]))
})

test_that("Gaussian subject null calibration is independent of the fitting code", {
  set.seed(661); s <- 20; r <- 2; p <- 1000
  ss <- group_fixture(s, r, p)
  for (j in seq_len(s)) ss[[j]]$estimate <- matrix(rnorm(p*r), p, r)
  g <- pattern_group(ss, ss[[1]]$basis)
  expect_lt(abs(mean(g$omnibus$p < .05) - .05), .025)
  expect_lt(abs(mean(g$omnibus$p) - .5), .035)
})

test_that("group prediction summaries retain unavailable subjects and reference identity", {
  ss <- group_fixture()
  for (j in seq_along(ss)) ss[[j]]$prediction <- data.frame(response = "y1", metric = "MSE", value = j)
  g <- pattern_group(ss, ss[[1]]$basis)
  expect_equal(g$prediction$summary$mean, 6.5)
  expect_equal(g$prediction$summary$n_subjects, 12)
  expect_output(print(g), "12 subjects, 3 features, 2 components")
  ss[[1]]$prediction$value <- NA_real_
  g <- pattern_group(ss, ss[[1]]$basis)
  expect_true(is.na(g$prediction$summary$mean))
  expect_equal(g$prediction$summary$n_subjects, 11)
  ref <- ss[[1]]$basis; ref$matrix <- ref$matrix %*% matrix(c(0, 1, -1, 0), 2)
  rotated <- pattern_group(ss, ref)
  expect_false(identical(rotated$reference_basis$basis_id, g$reference_basis$basis_id))
})

test_that("incompatible subspaces and corrupt uncertainty cannot be pooled", {
  ss <- group_fixture(); ref <- ss[[1]]$basis
  ref$matrix <- rbind(ref$matrix, 0); ref$target_ids <- c(ref$target_ids, "y3")
  for (j in seq_along(ss)) ss[[j]]$basis <- ref
  bad <- ss; bad[[1]]$basis$matrix[3, 1] <- 1
  expect_error(pattern_group(bad, ref), "same target subspace")
  bad <- ss; bad[[1]]$basis$type <- "categorical"
  expect_error(pattern_group(bad, ref), "incompatible")
  bad <- ss; bad[[1]]$covariance[1, 2, 1] <- 100
  expect_error(pattern_group(bad, ref), "symmetric positive semidefinite")
  bad <- ss; bad[[1]]$covariance[, , 1] <- 0
  expect_error(pattern_group(bad, ref, effects = "fixed"), "nonsingular")
  g <- pattern_group(bad, ref)
  expect_true(is.na(g$heterogeneity$Q[1]))
  mapping <- lapply(ss, function(x) setNames(x$feature_ids, x$feature_ids))
  names(mapping) <- vapply(ss, function(x) x$provenance$subject_id, character(1))
  g2 <- pattern_group(ss, ref, rev(mapping))
  expect_equal(g2$estimate, pattern_group(ss, ref)$estimate)
  names(mapping)[1] <- "unknown"
  expect_error(pattern_group(ss, ref, mapping), "match subject IDs")
  expect_error(pattern_group(ss, ref, list()), "one named character")
  expect_error(pattern_group(ss[1], ref), "at least two")
})

test_that("covariance validity is relative to measurement units", {
  ss <- group_fixture(); ref <- ss[[1]]$basis
  for (j in seq_along(ss)) {
    ss[[j]]$estimate <- ss[[j]]$estimate * 1e-10
    ss[[j]]$covariance <- ss[[j]]$covariance * 1e-20
  }
  expect_true(all(is.finite(pattern_group(ss, ref)$omnibus$p)))
  ss[[1]]$covariance[, , 1] <- diag(c(-1, 1)) * 1e-22
  expect_error(pattern_group(ss, ref), "positive semidefinite")
})

test_that("rank-three subject confirmation transports full sandwich covariance", {
  set.seed(667); n <- 80; p <- 8; r <- 3
  Y <- matrix(rnorm(n*r), n, r, dimnames = list(NULL, paste0("y", 1:r)))
  X <- matrix(rnorm(n*p), n, p)
  fit <- .pattern_fit(X, Y, rank = r, control = pattern_control())
  block <- rep(1:10, each = 8)
  test <- X + Y %*% matrix(rnorm(r*p), r, p)
  space <- neuroim2::NeuroSpace(c(2, 2, 2), c(1, 1, 1))
  dataset <- mvpa_dataset(neuroim2::NeuroVec(array(t(test), c(2, 2, 2, n)),
    neuroim2::NeuroSpace(c(2, 2, 2, n), c(1, 1, 1))),
    mask = neuroim2::NeuroVol(array(1, c(2, 2, 2)), space))
  design <- mvpa_design(data.frame(id = 1:n), cv_labels = factor(rep(1:2, n/2)), targets = Y)
  c <- pattern_confirm(fit, dataset, design, block,
    confirmation_plan("block_robust"), paste0("c", 1:n), paste0("d", 1:n),
    feature_ids = paste0("v", 1:p), preprocessing_id = "BOLD", subject_id = "s1")
  D <- cbind(.pattern_targets_apply(fit$y_transform, Y) %*% fit$C, 1)
  lmfit <- lm(test[, 1] ~ D - 1)
  bread <- solve(crossprod(D))
  meat <- crossprod(rowsum(D * residuals(lmfit), block))
  V <- bread %*% meat %*% bread * 10/9 * 79/76
  expect_equal(.pattern_confirmation_cov(c, 1), V[1:3, 1:3], ignore_attr = TRUE, tolerance = 1e-10)
  beta <- coef(lmfit)[1:3]
  expect_equal(c$omnibus$F[1], as.numeric(crossprod(beta, solve(V[1:3, 1:3], beta))/3), tolerance = 1e-10)
})
