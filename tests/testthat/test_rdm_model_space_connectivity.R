rdm_model_space_vec_to_mat <- function(v, n, labels = paste0("i", seq_len(n))) {
  M <- matrix(0, n, n, dimnames = list(labels, labels))
  M[lower.tri(M)] <- v
  M + t(M)
}

test_that("rdm_model_space_connectivity reduces to outer RSA scores for one model", {
  set.seed(101)
  p <- 30L
  a <- rnorm(p)
  roi <- cbind(
    roi1 = a + rnorm(p, sd = 0.2),
    roi2 = -0.7 * a + rnorm(p, sd = 0.4),
    roi3 = rnorm(p)
  )
  model <- matrix(a, ncol = 1L, dimnames = list(NULL, "A1"))

  got <- rdm_model_space_connectivity(roi, model, method = "pearson")

  r <- stats::cor(roi, a)
  expected <- tcrossprod(drop(r))
  dimnames(expected) <- list(colnames(roi), colnames(roi))

  expect_s3_class(got, "rdm_model_space_connectivity")
  expect_equal(got$rank, 1L)
  expect_equal(got$similarity, expected, tolerance = 1e-12)
  expect_equal(got$scores[, 1], drop(r), tolerance = 1e-12)
})

test_that("rdm_model_space_connectivity decomposes raw ROI similarity", {
  set.seed(102)
  p <- 45L
  latent <- matrix(rnorm(p * 4), p, 4)
  model <- latent %*% matrix(c(
    1, 0.6, 0.2, 0.1,
    0.4, 1, 0.3, 0.2,
    0.2, 0.4, 1, 0.5,
    0.1, 0.2, 0.5, 1
  ), 4, 4)
  colnames(model) <- paste0("A", 1:4)

  roi <- model %*% matrix(rnorm(4 * 6), 4, 6) +
    matrix(rnorm(p * 6, sd = 0.3), p, 6)
  colnames(roi) <- paste0("ROI", 1:6)

  got <- rdm_model_space_connectivity(roi, model, method = "pearson")

  component_sum <- Reduce(`+`, got$component_similarity)
  expect_equal(component_sum, got$similarity, tolerance = 1e-10)
  expect_equal(got$raw_similarity, got$similarity + got$residual_similarity, tolerance = 1e-10)
  expect_equal(unname(diag(got$profile_similarity)), rep(1, ncol(roi)), tolerance = 1e-10)
  expect_equal(rownames(got$model_axis_cor), colnames(got$scores))
  expect_equal(colnames(got$model_axis_cor), colnames(model))
})

test_that("rdm_model_space_connectivity accepts feature-RSA style vector tibbles", {
  set.seed(103)
  p <- 10L
  model <- matrix(rnorm(p * 2), p, 2, dimnames = list(NULL, c("common", "change")))
  roi_tbl <- tibble::tibble(
    roinum = c(10L, 20L, 30L),
    observation_index = list(1:5, 1:5, 1:5),
    rdm_vec = list(
      model[, 1] + rnorm(p, sd = 0.2),
      model[, 1] - model[, 2] + rnorm(p, sd = 0.2),
      -model[, 2] + rnorm(p, sd = 0.2)
    )
  )

  got <- rdm_model_space_connectivity(roi_tbl, model, method = "pearson")

  expect_equal(rownames(got$similarity), as.character(roi_tbl$roinum))
  expect_equal(colnames(got$similarity), as.character(roi_tbl$roinum))
  expect_equal(colnames(got$raw_model_scores), colnames(model))
})

test_that("rdm_model_space_connectivity accepts model RDM matrices and dist objects", {
  set.seed(104)
  n <- 6L
  p <- n * (n - 1L) / 2L
  model_vec <- matrix(rnorm(p * 2), p, 2)
  model_mats <- list(
    semantic = rdm_model_space_vec_to_mat(model_vec[, 1], n = n),
    visual = stats::as.dist(rdm_model_space_vec_to_mat(model_vec[, 2], n = n))
  )
  roi <- model_vec %*% matrix(c(1, 0.4, -0.2, 0.9), 2, 2) +
    matrix(rnorm(p * 2, sd = 0.2), p, 2)
  colnames(roi) <- c("left", "right")

  got <- rdm_model_space_connectivity(roi, model_mats)
  expected_scores <- stats::cor(roi, model_vec)
  dimnames(expected_scores) <- list(colnames(roi), names(model_mats))

  expect_equal(got$model_labels, names(model_mats))
  expect_equal(dim(got$similarity), c(2L, 2L))
  expect_equal(got$raw_model_scores, expected_scores, tolerance = 1e-12)
})

test_that("rdm_model_space_connectivity handles missing cells by complete rows", {
  set.seed(105)
  model <- matrix(rnorm(24), 12, 2, dimnames = list(NULL, c("A1", "A2")))
  roi <- matrix(rnorm(36), 12, 3, dimnames = list(NULL, paste0("R", 1:3)))
  roi[2, 1] <- NA_real_
  model[5, 2] <- NA_real_

  got <- rdm_model_space_connectivity(roi, model, use = "complete.obs")
  expect_equal(got$n_pairs_total, 12L)
  expect_equal(got$n_pairs_used, 10L)

  expect_error(
    rdm_model_space_connectivity(roi, model, use = "everything"),
    "finite"
  )
})
