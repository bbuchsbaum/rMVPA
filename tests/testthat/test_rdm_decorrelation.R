rdm_vec_to_mat <- function(v, n, labels = letters[seq_len(n)]) {
  M <- matrix(0, n, n, dimnames = list(labels, labels))
  M[lower.tri(M)] <- v
  M + t(M)
}

make_correlated_rdms <- function(n = 9, seed = 1) {
  set.seed(seed)
  m <- n * (n - 1L) / 2L
  z1 <- rnorm(m)
  z2 <- 0.75 * z1 + rnorm(m, sd = 0.5)
  z3 <- 0.55 * z2 + rnorm(m, sd = 0.6)
  z4 <- 0.45 * z3 + rnorm(m, sd = 0.7)
  vals <- list(z1, z2, z3, z4)
  names(vals) <- c("vgg4", "vgg7", "vgg11", "vgg16")
  lapply(vals, rdm_vec_to_mat, n = n)
}

test_that("identity shrinkage leaves rank RSA correlations unchanged", {
  set.seed(10)
  X1 <- matrix(rnorm(12 * 5), 12, 5)
  X2 <- matrix(rnorm(12 * 5), 12, 5)
  K1 <- stats::cor(t(X1))
  K2 <- stats::cor(t(X2))
  lambda <- 0.4
  K1_shrunk <- (1 - lambda) * K1 + lambda * diag(nrow(K1))
  K2_shrunk <- (1 - lambda) * K2 + lambda * diag(nrow(K2))

  d1 <- (1 - K1)[lower.tri(K1)]
  d2 <- (1 - K2)[lower.tri(K2)]
  s1 <- (1 - K1_shrunk)[lower.tri(K1_shrunk)]
  s2 <- (1 - K2_shrunk)[lower.tri(K2_shrunk)]

  expect_equal(
    stats::cor(d1, d2, method = "spearman"),
    stats::cor(s1, s2, method = "spearman"),
    tolerance = 1e-12
  )
  expect_equal(
    stats::cor(d1, d2, method = "kendall"),
    stats::cor(s1, s2, method = "kendall"),
    tolerance = 1e-12
  )
})

test_that("ordered innovation full residuals are orthogonal in RDM-vector space", {
  rdms <- make_correlated_rdms(n = 10, seed = 11)

  dec <- rdm_decorrelate(
    rdms,
    method = "ordered_innovation",
    similarity = "pearson",
    objective = "full_residual",
    return = "matrix"
  )

  C <- dec$cross_rdm_cor_after
  expect_lt(max(abs(C[lower.tri(C)])), 1e-10)
  expect_equal(unname(dec$gamma), c(1, 0, 0, 0))
  expect_true(all(vapply(dec$rdms, is.matrix, logical(1))))
  expect_true(all(vapply(dec$rdms, function(M) isTRUE(all.equal(M, t(M))), logical(1))))
})

test_that("constraint objective satisfies epsilon when feasible and preserves more when relaxed", {
  rdms <- make_correlated_rdms(n = 9, seed = 12)

  tight <- rdm_decorrelate(
    rdms,
    similarity = "pearson",
    epsilon = 0.02,
    gamma_grid = seq(0, 1, by = 0.1),
    return = "vector"
  )
  relaxed <- rdm_decorrelate(
    rdms,
    similarity = "pearson",
    epsilon = 0.5,
    gamma_grid = seq(0, 1, by = 0.1),
    return = "vector"
  )

  expect_lte(tight$mean_abs_cor_after, 0.02 + 1e-12)
  expect_gte(mean(relaxed$preservation), mean(tight$preservation) - 1e-12)
  expect_lte(relaxed$mean_abs_cor_after, 0.5 + 1e-12)
})

test_that("adjusted dist outputs plug into rsa_design", {
  rdms <- make_correlated_rdms(n = 8, seed = 120)
  dec <- rdm_decorrelate(
    rdms,
    similarity = "pearson",
    epsilon = 0.05,
    gamma_grid = seq(0, 1, by = 0.1)
  )

  des <- rsa_design(
    ~ vgg4 + vgg7 + vgg11 + vgg16,
    data = dec$rdms
  )

  expect_s3_class(dec$rdms$vgg4, "dist")
  expect_s3_class(des, "rsa_design")
  expect_equal(names(des$model_mat), names(rdms))
})

test_that("ZCA decorrelation whitens full-rank RDM vectors", {
  rdms <- make_correlated_rdms(n = 11, seed = 13)

  dec <- rdm_decorrelate(
    rdms,
    method = "zca",
    similarity = "pearson",
    return = "vector"
  )

  cov_after <- crossprod(as.matrix(dec$vectors)) / (nrow(dec$vectors) - 1L)
  expect_equal(unname(cov_after), diag(4), tolerance = 1e-8)
  expect_lt(dec$mean_abs_cor_after, 1e-8)
})

test_that("ordered decorrelation is invariant to common item permutation", {
  rdms <- make_correlated_rdms(n = 9, seed = 14)
  dec <- rdm_decorrelate(
    rdms,
    similarity = "pearson",
    objective = "full_residual",
    return = "matrix"
  )

  perm <- sample(seq_len(9))
  rdms_perm <- lapply(rdms, function(M) M[perm, perm])
  dec_perm <- rdm_decorrelate(
    rdms_perm,
    similarity = "pearson",
    objective = "full_residual",
    return = "matrix"
  )
  inv <- order(perm)
  unpermuted <- lapply(dec_perm$rdms, function(M) M[inv, inv])

  for (nm in names(rdms)) {
    expect_equal(unpermuted[[nm]], dec$rdms[[nm]], tolerance = 1e-10)
  }
  expect_equal(dec_perm$mean_abs_cor_after, dec$mean_abs_cor_after, tolerance = 1e-12)
})

test_that("correlation projection returns valid correlation-cone RDMs", {
  rdms <- make_correlated_rdms(n = 8, seed = 15)

  dec <- rdm_decorrelate(
    rdms,
    similarity = "pearson",
    objective = "full_residual",
    return = "matrix",
    project_correlation = TRUE
  )

  for (M in dec$rdms) {
    K <- 1 - M
    expect_equal(unname(diag(K)), rep(1, nrow(K)), tolerance = 1e-12)
    expect_true(isTRUE(all.equal(K, t(K), tolerance = 1e-10)))
    expect_gte(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  }
})

test_that("rdm_decorrelate validates malformed and degenerate RDM inputs", {
  rdms <- make_correlated_rdms(n = 6, seed = 16)

  bad_sym <- rdms
  bad_sym[[1]][1, 2] <- bad_sym[[1]][1, 2] + 1
  expect_error(rdm_decorrelate(bad_sym), "symmetric")

  bad_na <- rdms
  bad_na[[2]][2, 1] <- NA_real_
  bad_na[[2]][1, 2] <- NA_real_
  expect_error(rdm_decorrelate(bad_na), "finite")

  constant <- rdms
  constant[[3]] <- matrix(0, 6, 6)
  expect_error(rdm_decorrelate(constant), "non-constant")
})

test_that("duplicated RDMs yield explicit degenerate innovation diagnostics", {
  rdms <- make_correlated_rdms(n = 7, seed = 17)
  rdms[[2]] <- rdms[[1]]

  dec <- rdm_decorrelate(
    rdms,
    similarity = "pearson",
    objective = "full_residual",
    return = "vector"
  )

  expect_true(dec$diagnostics$degenerate_innovation[[2]])
  expect_true(all(is.finite(dec$preservation[-2])))
  expect_equal(dec$vectors[, 2], rep(0, nrow(dec$vectors)))
})

test_that("rdm_decorrelate representative problem stays within perf budget", {
  skip_on_cran()
  skip_if_not_perf_tests()

  n_items <- 32L
  rdms <- make_correlated_rdms(n = n_items, seed = 18)
  elapsed <- as.numeric(system.time({
    dec <- rdm_decorrelate(
      rdms,
      similarity = "spearman",
      epsilon = 0.05,
      gamma_grid = seq(0, 1, by = 0.05),
      return = "vector"
    )
  })["elapsed"])
  budget <- as.numeric(Sys.getenv("RMVPA_RDM_DECORRELATE_MAX_SECONDS", "5"))

  message(sprintf(
    "rdm_decorrelate perf guardrail: elapsed=%.3fs budget=%.3fs n=%d L=%d",
    elapsed, budget, n_items, length(rdms)
  ))
  expect_lte(elapsed, budget)
  expect_lte(dec$mean_abs_cor_after, 0.05 + 1e-12)
})
