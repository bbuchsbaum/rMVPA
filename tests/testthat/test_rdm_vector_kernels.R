test_that("shared RDM pair index cache matches lower.tri order", {
  info <- rMVPA:::.rdm_pair_indices(4L)
  expected <- which(lower.tri(matrix(TRUE, 4L, 4L)), arr.ind = TRUE)

  expect_equal(info$pair_idx, expected)
  expect_equal(info$i, expected[, 1])
  expect_equal(info$j, expected[, 2])
  expect_identical(rMVPA:::.crossnobis_pair_indices(4L), info)
})

test_that("correlation RDM vector kernels match full matrix references", {
  set.seed(8201)
  X <- matrix(rnorm(6 * 9), nrow = 6)
  rownames(X) <- paste0("C", seq_len(nrow(X)))

  ref_pearson <- rMVPA:::pairwise_dist(cordist(method = "pearson"), X)
  ref_spearman <- rMVPA:::pairwise_dist(cordist(method = "spearman"), X)

  expect_equal(
    rMVPA:::.rdm_vector_correlation(X, method = "pearson"),
    ref_pearson[lower.tri(ref_pearson)],
    tolerance = 1e-12
  )
  expect_equal(
    rMVPA:::.rdm_vector_correlation(X, method = "spearman"),
    ref_spearman[lower.tri(ref_spearman)],
    tolerance = 1e-12
  )
  expect_equal(
    rMVPA:::pairwise_dist_vector(cordist(method = "pearson"), X),
    ref_pearson[lower.tri(ref_pearson)],
    tolerance = 1e-12
  )
})

test_that("squared Euclidean RDM vector kernel matches public distance matrix", {
  set.seed(8202)
  X <- matrix(rnorm(7 * 5), nrow = 7)
  rownames(X) <- paste0("C", seq_len(nrow(X)))

  ref <- rMVPA:::pairwise_dist(eucdist(), X)
  expect_equal(
    rMVPA:::.rdm_vector_euclidean(X),
    ref[lower.tri(ref)],
    tolerance = 1e-12
  )
  expect_equal(
    rMVPA:::.rdm_vector_sqeuclidean(X),
    ref[lower.tri(ref)]^2,
    tolerance = 1e-12
  )
  expect_equal(
    rMVPA:::.rdm_vector_sqeuclidean(X, normalize_by_features = TRUE),
    ref[lower.tri(ref)]^2 / ncol(X),
    tolerance = 1e-12
  )
})

test_that("supplied-precision Mahalanobis RDM vector matches whitened Euclidean reference", {
  set.seed(8203)
  X <- matrix(rnorm(5 * 6), nrow = 5)
  A <- matrix(rnorm(ncol(X)^2), nrow = ncol(X))
  precision <- crossprod(A) + diag(ncol(X)) * 0.25
  W <- t(chol(precision))

  ref <- rMVPA:::.rdm_vector_sqeuclidean(X %*% W, normalize_by_features = TRUE)
  got <- rMVPA:::.rdm_vector_sqmahalanobis(
    X,
    precision = precision,
    normalize_by_features = TRUE
  )

  expect_equal(got, ref, tolerance = 1e-10)
  expect_equal(
    rMVPA:::pairwise_dist_vector(mahadist(), X, precision = precision),
    sqrt(rMVPA:::.rdm_vector_sqmahalanobis(X, precision = precision)),
    tolerance = 1e-10
  )
})

test_that("RDM vector kernels handle zero or one condition", {
  X0 <- matrix(numeric(0), nrow = 0, ncol = 3)
  X1 <- matrix(1:3, nrow = 1)

  expect_equal(rMVPA:::.rdm_vector_correlation(X0), numeric(0))
  expect_equal(rMVPA:::.rdm_vector_sqeuclidean(X0), numeric(0))
  expect_equal(rMVPA:::.rdm_vector_correlation(X1), numeric(0))
  expect_equal(rMVPA:::.rdm_vector_sqeuclidean(X1), numeric(0))
})
