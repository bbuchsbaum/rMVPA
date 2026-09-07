library(testthat)

# Dense reference of Psi from a noise object.
# Reference implementation of the sufficient statistic u = A' Psi^{-1} x.
.pattern_sufficient_ref <- function(fit, X) {
  Xc <- sweep(X[, fit$feature_index, drop = FALSE], 2L, fit$x_transform$mu, "-")
  Xc %*% fit$precision_A
}

.dense_psi <- function(noise) {
  Psi <- diag(noise$D)
  if (!is.null(noise$U)) Psi <- Psi + noise$U %*% t(noise$U)
  Psi
}

test_that("noise model identities hold against dense references", {
  set.seed(1)
  p <- 40; h <- 3
  D <- runif(p, 0.5, 2); U <- matrix(rnorm(p * h), p, h) * 0.7
  nz <- rMVPA:::new_pattern_noise(D, U)
  Psi <- .dense_psi(nz)
  Pinv <- solve(Psi)
  M <- matrix(rnorm(p * 5), p, 5)

  expect_equal(rMVPA:::.noise_apply_precision(nz, M), Pinv %*% M, tolerance = 1e-10)
  expect_equal(rMVPA:::.noise_apply(nz, M), Psi %*% M, tolerance = 1e-10)
  expect_equal(rMVPA:::.noise_precision_diag(nz), diag(Pinv), tolerance = 1e-10)
  # W' W = Psi^{-1} and W^{-1} inverts W
  W <- rMVPA:::.noise_whiten(nz, diag(p))
  expect_equal(crossprod(W), Pinv, tolerance = 1e-10)
  expect_equal(rMVPA:::.noise_unwhiten(nz, W), diag(p), tolerance = 1e-10)
  # quadratic form
  R <- matrix(rnorm(7 * p), 7, p)
  expect_equal(rMVPA:::.noise_quadform(nz, R), sum(diag(R %*% Pinv %*% t(R))), tolerance = 1e-9)
  # exact marginal restriction (not a submatrix of the precision)
  keep <- c(3, 9, 10, 25, 33)
  nr <- rMVPA:::.noise_restrict(nz, keep)
  expect_equal(rMVPA:::.noise_apply_precision(nr, diag(length(keep))), solve(Psi[keep, keep]), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(solve(Psi[keep, keep]), Pinv[keep, keep])))
  # diagonal-only object degrades gracefully
  nd <- rMVPA:::new_pattern_noise(D, NULL)
  expect_equal(rMVPA:::.noise_apply_precision(nd, M), M / D)
  expect_equal(nd$h, 0L)
})

test_that("noise estimation recovers a low-rank residual component and respects type", {
  set.seed(2)
  n <- 200; p <- 150
  U_true <- matrix(rnorm(p), p, 1) * 1.5
  E <- matrix(rnorm(n * p), n, p) + matrix(rnorm(n), n, 1) %*% t(U_true)
  nz <- rMVPA:::.estimate_pattern_noise(E, type = "diag_lowrank", rank = "auto", max_rank = 5)
  expect_equal(nz$h, 1L)
  # recovered direction aligns with the truth
  expect_gt(abs(cor(nz$U[, 1], U_true[, 1])), 0.95)
  # the total variance is preserved (D + diag(UU') ~ residual variance)
  expect_equal(nz$D + rowSums(nz$U^2), (1 - 0.1) * colSums(E^2) / (n - 1) + 0.1 * median(colSums(E^2) / (n - 1)),
               tolerance = 0.05)
  expect_equal(rMVPA:::.estimate_pattern_noise(E, type = "diag")$h, 0L)
  expect_equal(rMVPA:::.estimate_pattern_noise(E, type = "identity")$D, rep(1, p))
  expect_equal(rMVPA:::.estimate_pattern_noise(E, rank = 3)$h, 3L)
})

test_that("target coding whitens, caps rank at K - 1, and round-trips", {
  y <- factor(rep(c("a", "b", "c", "d"), 10))
  enc <- rMVPA:::.pattern_encode_targets(y)
  expect_equal(ncol(enc$Yw), 3L)
  expect_equal(crossprod(enc$Yw) / nrow(enc$Yw), diag(3), tolerance = 1e-10)
  expect_equal(enc$transform$priors, c(a = .25, b = .25, c = .25, d = .25))
  expect_equal(dim(enc$transform$class_codes), c(4L, 3L))
  # applying the transform to the same labels reproduces Yw
  expect_equal(rMVPA:::.pattern_targets_apply(enc$transform, y), enc$Yw)
  expect_error(rMVPA:::.pattern_targets_apply(enc$transform, factor("zzz")), "absent from training")

  set.seed(3)
  Y <- matrix(rnorm(50 * 3), 50, 3) %*% matrix(c(1, 2, 0, 0, 1, 1, 3, 0, 1), 3, 3) + 5
  encc <- rMVPA:::.pattern_encode_targets(Y, scale = "sd")
  expect_equal(crossprod(encc$Yw) / 50, diag(3), tolerance = 1e-10)
  expect_equal(unname(rMVPA:::.pattern_targets_invert(encc$transform, encc$Yw)), unname(Y), tolerance = 1e-10)
  # rank-deficient continuous targets keep only the non-null eigenspace
  Yd <- cbind(Y[, 1], Y[, 1] * 2, Y[, 2])
  encd <- rMVPA:::.pattern_encode_targets(Yd)
  expect_equal(encd$transform$q_eff, 2L)
})

test_that("unpenalized fit is the exact reduced-rank optimum and the alternation is monotone", {
  set.seed(4)
  n <- 80; p <- 60
  A <- matrix(0, p, 2); A[1:15, 1] <- rnorm(15); A[16:30, 2] <- rnorm(30 - 15)
  y <- factor(rep(c("a", "b", "c"), length.out = n))
  Yc <- scale(model.matrix(~ y - 1), scale = FALSE)
  C <- qr.Q(qr(matrix(rnorm(6), 3, 2)))
  X <- (Yc %*% C) %*% t(A) * 2 + matrix(rnorm(n * p), n, p)

  ctrl <- pattern_control(max_rank = 2, max_outer = 20)
  fit <- rMVPA:::.pattern_fit(X, y, rank = 2, control = ctrl)
  tr <- fit$diagnostics$objective
  expect_true(all(diff(tr) <= 1e-8 * abs(tr[1])))          # non-increasing
  expect_equal(tr[length(tr)], tr[1], tolerance = 1e-6)      # init already optimal
  expect_true(fit$diagnostics$converged)
  expect_equal(crossprod(fit$C), diag(2), tolerance = 1e-10) # orthonormal C
  expect_equal(fit$target_cov, diag(2), tolerance = 1e-8)    # Cov(T) = I under whitening
  # Gram identities
  expect_equal(fit$G, crossprod(fit$A, rMVPA:::.noise_apply_precision(fit$noise, fit$A)), tolerance = 1e-10)
  # recovered pattern subspace matches the truth (principal-angle cosines)
  cs <- svd(crossprod(qr.Q(qr(A)), qr.Q(qr(fit$A))))$d
  expect_true(all(cs > 0.9))
  # objective agrees with a dense evaluation
  Psi <- .dense_psi(fit$noise)
  Xc <- scale(X, scale = FALSE)
  Yw <- rMVPA:::.pattern_targets_apply(fit$y_transform, y)
  R <- Xc - (Yw %*% fit$C) %*% t(fit$A)
  expect_equal(tr[length(tr)], 0.5 * sum(diag(R %*% solve(Psi) %*% t(R))) / n, tolerance = 1e-8)
})

test_that("rank path is nested, rank is capped, and errors are informative", {
  set.seed(5)
  X <- matrix(rnorm(40 * 25), 40, 25)
  y <- factor(rep(c("a", "b", "c"), length.out = 40))
  path <- rMVPA:::.pattern_fit(X, y, rank = "path", control = pattern_control(max_rank = 5))
  expect_equal(names(path), c("rank1", "rank2"))               # K - 1 = 2
  expect_equal(path$rank1$A[, 1], path$rank2$A[, 1], tolerance = 1e-8)  # nested directions
  expect_error(rMVPA:::.pattern_fit(X, y, rank = 3), "exceeds the eligible rank")
  capped <- rMVPA:::.pattern_fit(X, y, rank = 3, cap_rank = TRUE)
  expect_equal(capped$rank, 2L)
  expect_error(rMVPA:::.pattern_fit(X[1:2, ], y[1:2], rank = 1), "at least three")
  expect_error(rMVPA:::.pattern_fit(X, factor(rep("a", 40)), rank = 1), "at least two classes")
})

test_that("constant and non-finite columns are dropped but positions are kept", {
  set.seed(6)
  X <- matrix(rnorm(30 * 12), 30, 12)
  X[, 4] <- 7; X[, 9] <- NA
  y <- factor(rep(c("a", "b"), 15))
  fit <- rMVPA:::.pattern_fit(X, y, rank = 1)
  expect_equal(fit$feature_index, setdiff(1:12, c(4, 9)))
  expect_equal(nrow(fit$A), 10L)
  expect_equal(fit$p_input, 12L)
  # prediction accepts either the full input width or the retained width
  p_full <- predict(fit, X, type = "prob")
  p_kept <- predict(fit, X[, fit$feature_index], type = "prob")
  expect_equal(p_full, p_kept)
  expect_error(predict(fit, X[, 1:5]), "expected 12")
})

test_that("prediction: probabilities, classes, decoding, encoding, scores, and degeneracy", {
  set.seed(7)
  n <- 90; p <- 50
  A <- matrix(0, p, 2); A[1:10, 1] <- 2; A[11:20, 2] <- -2
  y <- factor(rep(c("a", "b", "c"), 30))
  Yc <- scale(model.matrix(~ y - 1), scale = FALSE)
  C <- qr.Q(qr(matrix(rnorm(6), 3, 2)))
  X <- (Yc %*% C) %*% t(A) + matrix(rnorm(n * p), n, p)
  fit <- rMVPA:::.pattern_fit(X, y, rank = 2)

  P <- predict(fit, X, type = "prob")
  expect_equal(unname(rowSums(P)), rep(1, n))
  expect_equal(colnames(P), c("a", "b", "c"))
  cl <- predict(fit, X, type = "class")
  expect_s3_class(cl, "factor")
  expect_gt(mean(cl == y), 0.8)   # chance is 1/3
  Z <- predict(fit, X, type = "scores")
  expect_equal(dim(Z), c(n, 2L))
  # The calibration identity W'A = I (W = Psi^{-1} A G^{-1}) is what makes
  # z = t + W'eps an unbiased score; it holds exactly.
  W <- fit$precision_A %*% solve(fit$G)
  expect_equal(crossprod(W, fit$A), diag(2), tolerance = 1e-10)
  expect_equal(Z, .pattern_sufficient_ref(fit, X) %*% MASS::ginv(fit$G), tolerance = 1e-8)
  # and the leading score tracks the leading target dimension
  Tt <- rMVPA:::.pattern_targets_apply(fit$y_transform, y) %*% fit$C
  expect_gt(abs(cor(Z[, 1], Tt[, 1])), 0.9)
  expect_error(predict(fit, X, type = "decode") , NA)
  expect_error(predict(fit, X, type = "encode"), "targets")

  # degenerate fit (no retained signal) returns priors exactly
  fit0 <- fit
  fit0$A[] <- 0; fit0$precision_A[] <- 0; fit0$G[] <- 0
  P0 <- predict(fit0, X, type = "prob")
  expect_equal(unname(P0[1, ]), unname(fit$y_transform$priors))
  expect_equal(nrow(unique(P0)), 1L)

  # continuous: encode then decode a clean target recovers it
  set.seed(8)
  Y <- matrix(rnorm(n * 3), n, 3)
  A3 <- matrix(0, p, 3); A3[1:10, 1] <- 2; A3[11:20, 2] <- -2; A3[21:30, 3] <- 1.5
  Xr <- Y %*% matrix(rnorm(9), 3, 3) %*% t(A3) + matrix(rnorm(n * p), n, p) * 0.1
  fitr <- rMVPA:::.pattern_fit(Xr, Y, rank = 3)
  Xhat <- predict(fitr, Xr, type = "encode", targets = Y)
  expect_equal(attr(Xhat, "feature_index"), fitr$feature_index)
  dec <- predict(fitr, Xhat, type = "decode")
  expect_gt(min(diag(cor(dec, Y))), 0.99)
  expect_error(predict(fitr, Xr, type = "prob"), "categorical")
  # degenerate continuous fit returns the target means
  fitr0 <- fitr; fitr0$A[] <- 0; fitr0$precision_A[] <- 0; fitr0$G[] <- 0
  d0 <- predict(fitr0, Xr, type = "decode")
  expect_equal(unname(d0[5, ]), unname(colMeans(Y)))
})

test_that("a pattern_fit serializes without training data and predicts identically", {
  set.seed(9)
  X <- matrix(rnorm(60 * 30), 60, 30)
  y <- factor(rep(c("a", "b"), 30))
  fit <- rMVPA:::.pattern_fit(X, y, rank = 1)
  expect_null(fit$X); expect_null(fit$Yw)
  fit2 <- unserialize(serialize(fit, NULL))
  expect_identical(predict(fit2, X, type = "prob"), predict(fit, X, type = "prob"))
  expect_lt(object.size(fit), object.size(X) * 3)
  expect_output(print(fit), "pattern_fit: rank 1")
})

test_that("NA targets are rejected instead of becoming a phantom class", {
  set.seed(40)
  X <- matrix(rnorm(30 * 10), 30, 10)
  y <- factor(c(rep(c("a", "b"), 14), NA, NA))
  expect_error(rMVPA:::.pattern_fit(X, y, rank = 1), "missing values")
  Y <- matrix(rnorm(30 * 2), 30, 2); Y[3, 1] <- NA
  expect_error(rMVPA:::.pattern_fit(X, Y, rank = 1), "missing values")
})

test_that("the initialization SVD uses the small matrix and is exactly equivalent", {
  set.seed(41)
  n <- 60; p <- 400
  y <- factor(rep(c("a", "b", "c"), length.out = n))
  X <- matrix(rnorm(n * p), n, p)
  X[, 1:40] <- X[, 1:40] + (as.integer(y) - 2) * 1.5
  fit <- rMVPA:::.pattern_fit(X, y, rank = 2)
  Yw <- rMVPA:::.pattern_targets_apply(fit$y_transform, y)
  Xc <- scale(X, scale = FALSE)
  Xw <- t(rMVPA:::.noise_whiten(fit$noise, t(Xc)))
  B <- crossprod(Yw, Xw) / n
  ref <- svd(Yw %*% B, nu = 0, nv = 2)     # the n x p factorization it replaces
  new <- svd(B, nu = 0, nv = 2)
  expect_equal(abs(new$v[, 1:2]), abs(ref$v[, 1:2]), tolerance = 1e-10)
  expect_equal(sqrt(n) * new$d[1:2], ref$d[1:2], tolerance = 1e-10)
  expect_equal(fit$diagnostics$init_singular_values, ref$d[1:2], tolerance = 1e-8)
})

test_that("the low-rank noise Gram path matches a direct SVD when p >> n", {
  set.seed(42)
  n <- 40; p <- 300
  E <- matrix(rnorm(n * p), n, p) + matrix(rnorm(n), n, 1) %*% t(matrix(rnorm(p), p, 1) * 1.5)
  gram <- rMVPA:::.estimate_pattern_noise(E, rank = 2, max_rank = 4)      # p > 2n: Gram route
  # direct reference on the same standardized residuals
  v <- colSums(E^2) / (n - 1); med <- median(v[v > 0])
  D <- pmax(0.9 * v + 0.1 * med, 1e-8 * med)
  sv <- svd(sweep(E, 2L, sqrt(D), "/"), nu = 0, nv = 4)
  expect_equal(abs(gram$U[, 1] / sqrt(D)) / sqrt(sum((gram$U[, 1] / sqrt(D))^2)),
               abs(sv$v[, 1]), tolerance = 1e-8)
  expect_equal(gram$meta$spectrum[1:2], sv$d[1:2]^2 / (n - 1), tolerance = 1e-8)
  expect_true(all(gram$D > 0))
})

test_that("pattern_control validates and normalizes its inputs", {
  ctrl <- pattern_control(max_rank = 3, noise = list(type = "diag"))
  expect_s3_class(ctrl, "pattern_control")
  expect_equal(ctrl$noise$type, "diag")
  expect_equal(ctrl$noise$max_rank, 10L)
  expect_error(pattern_control(max_rank = 0), "positive integer")
  expect_error(pattern_control(lambda_2 = -1), "must be 0")
  expect_error(pattern_control(lambda_2 = 0.1), "not implemented")
  expect_error(pattern_control(noise = list(type = "banana")))
})
