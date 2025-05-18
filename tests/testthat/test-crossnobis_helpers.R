context("compute_crossnobis_distances_sl order")

library(testthat)
library(rMVPA)

# Simple test to verify pair order mirrors lower.tri

test_that("compute_crossnobis_distances_sl returns pairs in lower.tri order", {
  K <- 3; V <- 1; M <- 2
  # create U_folds array filled with zeros
  U_folds <- array(0, dim = c(K, V, M),
                   dimnames = list(paste0("Cond", 1:K), NULL, NULL))
  res <- rMVPA:::compute_crossnobis_distances_sl(U_folds, P_voxels = V)

  pair_idx <- which(lower.tri(matrix(1, K, K)), arr.ind = TRUE)
  expected_names <- paste0(paste0("Cond", pair_idx[,1]), "_vs_", paste0("Cond", pair_idx[,2]))
  expect_equal(names(res), expected_names)
  expect_equal(length(res), K * (K - 1) / 2)
})
