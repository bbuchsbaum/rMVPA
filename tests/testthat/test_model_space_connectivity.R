test_that("model_space_connectivity.default forwards to rdm_model_space_connectivity", {
  set.seed(1)
  p <- 30L
  model <- matrix(rnorm(p * 2), p, 2, dimnames = list(NULL, c("M1", "M2")))
  roi   <- model %*% matrix(rnorm(2 * 4), 2, 4) +
           matrix(rnorm(p * 4, sd = 0.3), p, 4)
  colnames(roi) <- paste0("R", seq_len(4))

  expected <- rdm_model_space_connectivity(roi, model, method = "pearson")
  got      <- model_space_connectivity(roi, model, method = "pearson")

  expect_equal(got$similarity, expected$similarity, tolerance = 1e-12)
  expect_equal(got$scores,     expected$scores,     tolerance = 1e-12)
})

test_that("model_space_connectivity dispatches on regional_mvpa_result fingerprints", {
  set.seed(2)
  axis_names <- c("PC1", "PC2", "PC3")
  scores <- matrix(rnorm(15), 5, 3, dimnames = list(paste0("ROI", seq_len(5)), axis_names))
  fake <- structure(list(model_spec = NULL,
                         performance_table = tibble::tibble(roinum = seq_len(5)),
                         prediction_table = NULL,
                         vol_results = list(),
                         fits = NULL),
                    class = c("regional_mvpa_result", "list"))
  attr(fake, "fingerprints") <- list(ids = seq_len(5), scores = scores, rdms = NULL)

  conn <- model_space_connectivity(fake)
  expect_s3_class(conn, "rdm_model_space_connectivity")
  expect_equal(conn$scores, scores, tolerance = 1e-12)
  expect_equal(conn$similarity, tcrossprod(scores), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(rownames(conn$similarity), rownames(scores))
  expect_null(conn$raw_similarity)
  expect_null(conn$residual_similarity)
  expect_false(conn$decomposition_available)

  component_sum <- Reduce(`+`, conn$component_similarity)
  expect_equal(component_sum, conn$similarity, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(unname(diag(conn$profile_similarity)), rep(1, nrow(scores)),
               tolerance = 1e-12)
  expect_true(all(abs(conn$profile_similarity) <= 1 + 1e-12))
})

test_that("model_space_connectivity on a fingerprint-less searchlight result errors clearly", {
  fake <- structure(list(results = list(), metrics = character()),
                    class = c("searchlight_result", "list"))
  expect_error(model_space_connectivity(fake), "no stored fingerprints")
})

test_that("run_regional captures fingerprints end-to-end and feeds model_space_connectivity", {
  skip_if_not_installed("neuroim2")
  set.seed(404)
  n <- 16L
  arr  <- array(rnorm(prod(c(3, 3, 3, n))), c(3, 3, 3, n))
  sp   <- neuroim2::NeuroSpace(c(3, 3, 3, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(3, 3, 3)),
                                    neuroim2::NeuroSpace(c(3, 3, 3)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  R1 <- as.matrix(stats::dist(matrix(rnorm(n * 4), n, 4)))
  R2 <- as.matrix(stats::dist(matrix(rnorm(n * 4), n, 4)))
  des <- pair_rsa_design(items_a = paste0("it", seq_len(n)),
                         model = list(R1 = R1, R2 = R2))
  mdl <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson",
                   return_fingerprint = TRUE)

  region_idx <- array(rep(c(1L, 2L), each = 14L)[seq_len(27L)], c(3, 3, 3))
  region <- neuroim2::ClusteredNeuroVol(mask,
                                        as.integer(region_idx[mask > 0]))

  res <- run_regional(mdl, region, verbose = FALSE)
  fp  <- attr(res, "fingerprints")
  expect_false(is.null(fp))
  expect_equal(dim(fp$scores), c(2L, 2L))
  expect_setequal(fp$ids, c(1L, 2L))

  conn <- model_space_connectivity(res)
  expect_s3_class(conn, "rdm_model_space_connectivity")
  expect_equal(dim(conn$similarity), c(2L, 2L))
  expect_equal(conn$similarity, tcrossprod(fp$scores), tolerance = 1e-12,
               ignore_attr = TRUE)
})

test_that("fingerprint-only connectivity is permutation equivariant", {
  set.seed(44)
  scores <- matrix(rnorm(24), 6, 4,
                   dimnames = list(paste0("ROI", 1:6), paste0("PC", 1:4)))
  perm <- c(3L, 6L, 1L, 5L, 2L, 4L)

  make_result <- function(scores, ids) {
    fake <- structure(list(model_spec = NULL,
                           performance_table = tibble::tibble(roinum = ids),
                           prediction_table = NULL,
                           vol_results = list(),
                           fits = NULL),
                      class = c("regional_mvpa_result", "list"))
    attr(fake, "fingerprints") <- list(ids = ids, scores = scores, rdms = NULL)
    fake
  }

  conn <- model_space_connectivity(make_result(scores, seq_len(6)))
  conn_perm <- model_space_connectivity(make_result(scores[perm, , drop = FALSE], perm))

  expect_equal(conn_perm$similarity[rownames(conn$similarity),
                                    colnames(conn$similarity)],
               conn$similarity, tolerance = 1e-12)
  expect_equal(conn_perm$profile_similarity[rownames(conn$profile_similarity),
                                            colnames(conn$profile_similarity)],
               conn$profile_similarity, tolerance = 1e-12)
})

test_that("rsa_model fingerprint matches rdm_model_space_connectivity scores", {
  skip_if_not_installed("neuroim2")
  set.seed(3)
  n <- 14L
  v <- 25L
  X <- matrix(rnorm(n * v), n, v)
  R1 <- as.matrix(stats::dist(matrix(rnorm(n * 3), n, 3)))
  R2 <- as.matrix(stats::dist(matrix(rnorm(n * 3), n, 3)))

  des <- rsa_design(~ R1 + R2, list(R1 = R1, R2 = R2))
  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  m <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson",
                 return_fingerprint = TRUE,
                 fingerprint_method = "pearson", fingerprint_basis = "pca")

  out <- train_model(m, X, y = NULL, indices = NULL)
  fp  <- attr(out, "fingerprint", exact = TRUE)
  expect_true(is.numeric(fp))
  expect_equal(length(fp), 2L)  # two model RDMs => rank 2

  # Reference: compute the neural pair vector, then run rdm_model_space_connectivity
  dtrain <- 1 - stats::cor(t(X), method = "pearson")
  dvec   <- dtrain[lower.tri(dtrain)]

  model_mat <- do.call(cbind, des$model_mat)
  ref <- rdm_model_space_connectivity(matrix(dvec, ncol = 1L,
                                              dimnames = list(NULL, "roi")),
                                       model_mat,
                                       method = "pearson", basis = "pca")
  expect_equal(unname(fp), unname(ref$scores[1, ]), tolerance = 1e-10)
})

test_that("rsa_model fingerprints exclude nuisance predictors", {
  skip_if_not_installed("neuroim2")
  set.seed(31)
  n <- 12L
  v <- 20L
  X <- matrix(rnorm(n * v), n, v)
  R <- as.matrix(stats::dist(matrix(rnorm(n * 3), n, 3)))
  N <- as.matrix(stats::dist(seq_len(n)))

  des <- pair_rsa_design(items_a = paste0("it", seq_len(n)),
                         model = list(signal = R),
                         nuisance = list(order = N))

  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  m <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson",
                 return_fingerprint = TRUE)
  out <- train_model(m, X, y = NULL, indices = NULL)
  fp <- attr(out, "fingerprint", exact = TRUE)

  expect_named(out, c("signal", "order"))
  expect_equal(length(fp), 1L)

  dtrain <- 1 - stats::cor(t(X), method = "pearson")
  dvec <- dtrain[lower.tri(dtrain)]
  ref <- rdm_model_space_connectivity(
    matrix(dvec, ncol = 1L, dimnames = list(NULL, "roi")),
    matrix(des$model_mat$signal, ncol = 1L, dimnames = list(NULL, "signal")),
    method = "pearson", basis = "pca"
  )
  expect_equal(unname(fp), unname(ref$scores[1, ]), tolerance = 1e-10)
})

test_that("classic rsa_design nuisance predictors are excluded from fingerprints", {
  skip_if_not_installed("neuroim2")
  set.seed(32)
  n <- 12L
  v <- 20L
  X <- matrix(rnorm(n * v), n, v)
  R <- as.matrix(stats::dist(matrix(rnorm(n * 3), n, 3)))
  N <- as.matrix(stats::dist(seq_len(n)))

  des <- rsa_design(
    ~ signal,
    data = list(signal = R),
    nuisance = list(order = N)
  )

  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  m <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson",
                 return_fingerprint = TRUE)
  out <- train_model(m, X, y = NULL, indices = NULL)
  fp <- attr(out, "fingerprint", exact = TRUE)

  expect_named(out, c("signal", "order"))
  expect_equal(length(fp), 1L)
  expect_identical(m$design$model_predictors, "signal")
  expect_identical(m$design$nuisance_predictors, "order")
})

test_that("rsa_model fingerprint basis drops near-collinear model axes without non-finite scores", {
  skip_if_not_installed("neuroim2")
  set.seed(515)
  n <- 10L
  X <- matrix(rnorm(n * 16), n, 16)
  R1 <- as.matrix(stats::dist(matrix(rnorm(n * 3), n, 3)))
  R2 <- R1

  arr  <- array(rnorm(prod(c(2, 2, 2, n))), c(2, 2, 2, n))
  sp   <- neuroim2::NeuroSpace(c(2, 2, 2, n))
  vec  <- neuroim2::NeuroVec(arr, sp)
  mask <- neuroim2::LogicalNeuroVol(array(1, c(2, 2, 2)),
                                    neuroim2::NeuroSpace(c(2, 2, 2)))
  ds   <- mvpa_dataset(train_data = vec, mask = mask)

  des <- pair_rsa_design(items_a = paste0("it", seq_len(n)),
                         model = list(R1 = R1, R2 = R2))
  m <- rsa_model(ds, des, distmethod = "pearson", regtype = "pearson",
                 return_fingerprint = TRUE)
  fp <- attr(train_model(m, X, y = NULL, indices = NULL), "fingerprint",
             exact = TRUE)

  expect_equal(length(fp), 1L)
  expect_true(all(is.finite(fp)))
})
