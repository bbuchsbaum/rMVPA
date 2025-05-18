library(testthat)
library(rMVPA)

context("compute_crossnobis_distances_sl")

# Test unbiasedness on simulated data

test_that("compute_crossnobis_distances_sl is unbiased on simple simulated data", {
  set.seed(123)
  V <- 5; K <- 2; M <- 3
  true_delta <- c(1, -0.5, 0.5, 2, -1)
  noise_sd <- 0.5
  nrep <- 200
  crossnobis_vals <- numeric(nrep)
  naive_vals <- numeric(nrep)

  for (r in seq_len(nrep)) {
    U_folds <- array(0, dim = c(K, V, M),
                     dimnames = list(paste0("C", 1:K), paste0("V", 1:V), NULL))
    for (m in seq_len(M)) {
      noise1 <- rnorm(V, 0, noise_sd)
      noise2 <- rnorm(V, 0, noise_sd)
      U_folds[1, , m] <- true_delta / 2 + noise1
      U_folds[2, , m] <- -true_delta / 2 + noise2
    }

    crossnobis_vals[r] <- compute_crossnobis_distances_sl(U_folds, P_voxels = V)[1]

    mean_patterns <- apply(U_folds, c(1, 2), mean)
    delta_mean <- mean_patterns[1, ] - mean_patterns[2, ]
    naive_vals[r] <- sum(delta_mean^2) / V
  }

  expected <- sum(true_delta^2) / V
  expect_lt(abs(mean(crossnobis_vals) - expected), 0.1)
  expect_gt(mean(naive_vals), expected)
})

# Test that pair ordering matches lower.tri

test_that("compute_crossnobis_distances_sl pair ordering matches lower.tri", {
  K <- 3; V <- 2; M <- 2
  U_folds <- array(0, dim = c(K, V, M),
                   dimnames = list(paste0("Cond", 1:K), NULL, NULL))
  U_folds[1, , ] <- 1
  U_folds[2, , ] <- 2
  U_folds[3, , ] <- 3

  dvec <- compute_crossnobis_distances_sl(U_folds, P_voxels = V)
  expect_equal(names(dvec), c("Cond2_vs_Cond1", "Cond3_vs_Cond1", "Cond3_vs_Cond2"))
})

