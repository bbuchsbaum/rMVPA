context("Joint off-diagonal endpoint decomposition")

.endpoint_fixture <- function(n = 12L, p = 15L, seed = 92L) {
  set.seed(seed)
  keys <- sprintf("item%02d", seq_len(n))
  E <- matrix(rnorm(n * p), n, dimnames = list(keys, NULL))
  R <- E + matrix(rnorm(n * p), n)
  rownames(R) <- keys
  list(E = E, R = R, keys = keys, block = rep("pool", n))
}

.endpoint_from_matrix <- function(M, block = rep("pool", nrow(M))) {
  era_endpoint_decomposition(similarity_matrix = M,
    item_key = paste0("item", seq_len(nrow(M))), competitor_block = block)
}

test_that("complete-block solution agrees with independent constrained lm", {
  set.seed(193)
  for (n in c(3L, 4L, 9L, 20L)) {
    M <- matrix(rnorm(n * n), n)
    z <- .endpoint_from_matrix(M)$items
    ij <- which(row(M) != col(M), arr.ind = TRUE)
    dat <- data.frame(y = M[ij], enc = factor(ij[, 1]), ret = factor(ij[, 2]))
    fit <- lm(y ~ enc + ret, dat, contrasts = list(enc = contr.sum, ret = contr.sum))
    cf <- coef(fit)
    a <- c(cf[seq.int(2L, n)], -sum(cf[seq.int(2L, n)]))
    b <- c(cf[seq.int(n + 1L, 2L * n - 1L)], -sum(cf[seq.int(n + 1L, 2L * n - 1L)]))
    expect_equal(z$shared, rep(unname(cf[1]), n), tolerance = 1e-12)
    expect_equal(z$encoding, unname(a), tolerance = 1e-12)
    expect_equal(z$retrieval, unname(b), tolerance = 1e-12)
    residual <- M - z$shared[1] - outer(z$encoding, z$retrieval, "+")
    diag(residual) <- 0
    expect_equal(rowSums(residual), rep(0, n), tolerance = 1e-12)
    expect_equal(colSums(residual), rep(0, n), tolerance = 1e-12)
    expect_equal(z$total, z$shared + z$encoding + z$retrieval + z$excess, tolerance = 1e-12)
  }
})

test_that("diagonal perturbation, transpose and specificity identities hold", {
  set.seed(17)
  n <- 7L
  a <- rnorm(n); a <- a - mean(a)
  b <- rnorm(n); b <- b - mean(b)
  M <- 0.2 + outer(a, b, "+")
  z <- .endpoint_from_matrix(M)$items
  expect_equal(z$excess, rep(0, n), tolerance = 1e-12)
  expect_equal(z$encoding, a, tolerance = 1e-12)
  expect_equal(z$retrieval, b, tolerance = 1e-12)
  injection <- rnorm(n)
  diag(M) <- diag(M) + injection
  z2 <- .endpoint_from_matrix(M)$items
  expect_equal(z2[, c("shared", "encoding", "retrieval")], z[, c("shared", "encoding", "retrieval")])
  expect_equal(z2$excess - z$excess, injection, tolerance = 1e-12)
  zt <- .endpoint_from_matrix(t(M))$items
  expect_equal(zt$encoding, z2$retrieval, tolerance = 1e-12)
  expect_equal(zt$retrieval, z2$encoding, tolerance = 1e-12)
  expect_equal(zt$excess, z2$excess, tolerance = 1e-12)
  off <- M; diag(off) <- 0
  expect_equal(diag(M) - colSums(off) / (n - 1), n / (n - 1) * z2$encoding + z2$excess)
  expect_equal(diag(M) - rowSums(off) / (n - 1) - colSums(off) / (n - 1) + z2$shared,
               z2$excess + (z2$encoding + z2$retrieval) / (n - 1))
})

test_that("fast pattern sums, explicit matrices and reordered keys agree", {
  for (seed in 1:8) {
    fx <- .endpoint_fixture(n = 17L, seed = seed)
    block <- setNames(rep(c("a", "b", "c"), c(3, 5, 9)), fx$keys)
    fast <- era_endpoint_decomposition(fx$E, fx$R, competitor_block = block)
    dense <- era_endpoint_decomposition(fx$E, fx$R, competitor_block = block, return_matrices = TRUE)
    oracle <- era_endpoint_decomposition(similarity_matrix = cor(t(fx$E), t(fx$R)), competitor_block = block)
    expect_equal(fast$items[, 3:7], dense$items[, 3:7], tolerance = 1e-12)
    expect_equal(fast$items[, 3:7], oracle$items[, 3:7], tolerance = 1e-12)
    expect_null(fast$matrices)
    expect_equal(length(dense$matrices), 3L)
    perm1 <- sample(17); perm2 <- sample(17)
    reordered <- era_endpoint_decomposition(fx$E[perm1, ], fx$R[perm2, ], competitor_block = block)
    expect_equal(unname(as.matrix(reordered$items[match(fx$keys, reordered$items$item_key), 3:7])),
                 unname(as.matrix(fast$items[, 3:7])), tolerance = 1e-12)
    explicit <- era_endpoint_decomposition(unname(fx$E[perm1, ]), unname(fx$R[perm2, ]),
      item_key = list(encoding = fx$keys[perm1], retrieval = fx$keys[perm2]), competitor_block = block)
    expect_equal(explicit$items, reordered$items)
    raw <- summarize_era_endpoints(fast)
    expect_equal(raw$blocks$raw_weight, c(3, 5, 9) / 17)
    expect_lt(raw$identity_error, 1e-12)
    expect_equal(raw$table$estimate[3:4], c(0, 0), tolerance = 1e-12)
    expect_true(all(raw$table$structural_zero[3:4]))
    expect_true(all(is.na(raw$table$p_value)))
  }
})

test_that("constant retrieval and encoding controls do not invent excess correlations", {
  set.seed(10)
  n <- 40L; p <- 80L
  vividness <- seq(-2, 2, length.out = n)
  q <- rnorm(p); q <- q - mean(q); q <- q / sqrt(sum(q^2))
  noise <- matrix(rnorm(n * p), n, p)
  noise <- noise - rowMeans(noise)
  noise <- noise - as.vector(noise %*% q) %o% q
  noise <- noise / sqrt(rowSums(noise^2))
  w <- 0.6 * tanh(vividness)
  E <- w %o% q + sqrt(1 - w^2) * noise
  R <- matrix(q, n, p, byrow = TRUE)
  for (transpose in c(FALSE, TRUE)) {
    z <- era_endpoint_decomposition(if (transpose) R else E, if (transpose) E else R,
      item_key = paste0("i", 1:n), competitor_block = rep("pool", n))
    expect_lt(max(abs(z$items$excess)), 1e-12)
    for (method in c("pearson", "spearman")) {
      assoc <- summarize_era_endpoints(z, "covariance_contributions", behavior = vividness, cor_method = method)
      tab <- assoc$table
      expect_true(is.na(tab$correlation[tab$quantity == "excess"]))
      expect_identical(tab$correlation_status[tab$quantity == "excess"], "constant_component")
      expect_gt(tab$contribution[1], 0.98)
      endpoint <- if (transpose) "retrieval" else "encoding"
      expect_equal(tab$contribution[tab$quantity == endpoint], tab$contribution[1], tolerance = 1e-12)
      expect_lt(assoc$identity_error, 1e-12)
    }
  }
})

test_that("common residualization and common cases preserve covariance additivity", {
  fx <- .endpoint_fixture(n = 24L)
  z <- era_endpoint_decomposition(fx$E, fx$R, competitor_block = rep(1:3, each = 8))
  set.seed(99)
  v <- rnorm(24); v[2] <- NA_real_
  motion <- matrix(rnorm(48), 24, 2); motion[4, 1] <- Inf
  for (adjust in c(FALSE, TRUE)) {
    out <- summarize_era_endpoints(z, "covariance_contributions", behavior = v,
                                  nuisance = motion, adjust_blocks = adjust)
    expect_equal(out$n, 22L)
    expect_equal(out$exclusions$reason[c(2, 4)], c("nonfinite_behavior", "nonfinite_nuisance"))
    keep <- out$exclusions$included
    X <- cbind(1, motion[keep, ])
    if (adjust) X <- cbind(X, model.matrix(~ factor(z$items$competitor_block[keep]))[, -1])
    d_res <- lm.fit(X, z$items$total[keep])$residuals
    v_res <- lm.fit(X, v[keep])$residuals
    expect_equal(out$table$contribution[1], cor(d_res, v_res), tolerance = 1e-12)
    expect_equal(sum(out$table$contribution[-1]), cor(d_res, v_res), tolerance = 1e-12)
    expect_lt(out$identity_error, 1e-12)
    expect_identical(out$shared_absorbed, adjust)
    if (!adjust) expect_gt(abs(out$table$contribution[2]), 1e-6)
    expect_gt(max(abs(out$subset_component_means[c("encoding", "retrieval")])), 1e-6)
    ranks <- summarize_era_endpoints(z, "covariance_contributions", behavior = v,
      nuisance = motion, adjust_blocks = adjust, cor_method = "spearman")
    expect_equal(ranks$table$contribution, out$table$contribution)
    expect_false(isTRUE(all.equal(ranks$table$correlation, out$table$correlation)))
    ranked_d <- lm.fit(X, rank(z$items$total[keep]))$residuals
    ranked_v <- lm.fit(X, rank(v[keep]))$residuals
    expect_equal(ranks$table$correlation[1], cor(ranked_d, ranked_v), tolerance = 1e-12)
    # Explicit intercept plus a redundant column spans the same nuisance space.
    redundant <- summarize_era_endpoints(z, "covariance_contributions", behavior = v,
      nuisance = cbind(1, motion, motion[, 1]), adjust_blocks = adjust)
    expect_equal(redundant$table$contribution, out$table$contribution, tolerance = 1e-12)
  }
  constant <- summarize_era_endpoints(z, "covariance_contributions", behavior = rep(1, 24))
  expect_identical(constant$status, "constant_behavior")
  expect_true(all(is.na(constant$table$contribution)))
  expect_true(all(is.na(constant$table$correlation)))
  saturated <- summarize_era_endpoints(z, "covariance_contributions", behavior = seq_len(24), nuisance = diag(24))
  expect_identical(saturated$status, "insufficient_residual_df")
  expect_true(all(is.na(saturated$table$contribution)))
})

test_that("unsupported geometry returns explicit block statuses or refuses", {
  fx <- .endpoint_fixture(n = 8L)
  block <- rep(c("small", "valid"), c(2, 6))
  z <- era_endpoint_decomposition(fx$E, fx$R, competitor_block = block)
  expect_equal(z$blocks$status, c("too_few_items", "ok"))
  expect_true(all(is.na(z$items$excess[1:2])))
  expect_equal(summarize_era_endpoints(z)$n, 6)
  expect_error(era_endpoint_decomposition(fx$E, fx$R, competitor_block = block, invalid = "error"), "too_few_items")
  for (kind in c("missing", "constant")) {
    E <- fx$E
    if (kind == "missing") E[3, 1] <- NA else E[3, ] <- 4
    bad <- era_endpoint_decomposition(E, fx$R, competitor_block = block)
    expect_true(all(is.na(bad$items$excess)))
    expect_identical(bad$blocks$status[2], if (kind == "missing") "nonfinite_patterns" else "constant_pattern")
    expect_equal(summarize_era_endpoints(bad)$n, 0)
  }
  scores <- matrix(rnorm(64), 8); scores[3, 4] <- NA
  expect_identical(.endpoint_from_matrix(scores, block)$blocks$status[2], "nonfinite_scores")
  scores[3, 4] <- 1; scores[1:2, 3:8] <- NA
  expect_identical(.endpoint_from_matrix(scores, block)$blocks$status[2], "ok")
  expect_error(era_endpoint_decomposition(fx$E, fx$R, competitor_block = block, pair_weights = matrix(1,8,8)), "Weighted")
  expect_error(era_endpoint_decomposition(unname(fx$E), unname(fx$R), competitor_block = block), "explicit item keys")
  expect_error(era_endpoint_decomposition(fx$E, fx$R, item_key = rep("a", 8), competitor_block = block), "unique")
  expect_error(era_endpoint_decomposition(fx$E, fx$R[-1, ], competitor_block = block), "one-to-one")
  expect_error(era_endpoint_decomposition(fx$E, fx$R, competitor_block = "all"), "one value per item")
  R <- fx$R; colnames(R) <- letters[1:ncol(R)]
  expect_error(era_endpoint_decomposition(fx$E, R, competitor_block = block), "feature columns")
  expect_error(.endpoint_from_matrix(matrix(1, 4, 3)), "square")
  z <- era_endpoint_decomposition(fx$E, fx$R, competitor_block = fx$block)
  expect_error(summarize_era_endpoints(z, behavior = seq_len(8)), "without behavior")
  expect_error(summarize_era_endpoints(z, "covariance_contributions", behavior = seq_len(7)), "one value per item")
})

test_that("small shared variation is retained in the additive covariance identity", {
  n <- 12L
  M <- matrix(0, n, n)
  M[1:6, 1:6] <- 0.1
  M[7:12, 7:12] <- 0.1 + 1e-8
  diag(M) <- diag(M) + seq_len(n) * 2e-8
  z <- .endpoint_from_matrix(M, rep(1:2, each = 6))
  out <- summarize_era_endpoints(z, "covariance_contributions", behavior = seq_len(n), adjust_blocks = FALSE)
  expect_gt(abs(out$table$contribution[2]), 0.01)
  expect_false(out$shared_absorbed)
  expect_lt(out$identity_error, 1e-8)
})

test_that("joint-excess legacy associations guard the actual complete-case subset", {
  keys <- letters[1:8]
  score <- c(0.1 + seq_len(4) * 1e-16, 0.2, 0.4, 0.1, 0.7)
  data <- data.frame(.era_key = keys, v = c(1:4, rep(NA_real_, 4)), row.names = keys)
  spec <- rMVPA:::.era_prepare_correlates(~ v, data)
  guarded <- rMVPA:::.era_correlate_metrics(score, keys, spec, "pearson", 3L,
                                           tolerance = sqrt(.Machine$double.eps))
  expect_true(is.na(guarded[["era_v_cor"]]))
  expect_equal(guarded[["era_v_n"]], 4L)
  reg <- rMVPA:::.era_prepare_association(~ v, ~ v, data)
  guarded_reg <- rMVPA:::.era_association_metrics(score, keys, reg, 3L,
                                                 tolerance = sqrt(.Machine$double.eps))
  expect_true(is.na(guarded_reg[["era_assoc_part_r_v"]]))
})
