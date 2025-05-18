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

# Helper to manually compute crossnobis distances
manual_crossnobis <- function(U_folds) {
  K <- dim(U_folds)[1]; V <- dim(U_folds)[2]; M <- dim(U_folds)[3]
  cond_pairs <- utils::combn(K, 2)
  n_pairs <- ncol(cond_pairs)
  res <- numeric(n_pairs)
  pair_names <- character(n_pairs)
  for (p in seq_len(n_pairs)) {
    i <- cond_pairs[1, p]; j <- cond_pairs[2, p]
    deltas <- U_folds[i, , ] - U_folds[j, , ]
    ip_mat <- crossprod(deltas)
    off_diag <- sum(ip_mat) - sum(diag(ip_mat))
    res[p] <- off_diag / (V * M * (M - 1))
    pair_names[p] <- paste0(dimnames(U_folds)[[1]][i], "_vs_", dimnames(U_folds)[[1]][j])
  }
  names(res) <- pair_names
  res
}


test_that("compute_crossnobis_distances_sl matches manual computation and is invariant to fold order", {
  set.seed(123)
  K <- 3; V <- 4; M <- 3
  U_folds <- array(rnorm(K * V * M), dim = c(K, V, M),
                   dimnames = list(paste0("C", 1:K), NULL, paste0("F", 1:M)))

  expected <- manual_crossnobis(U_folds)
  res <- compute_crossnobis_distances_sl(U_folds, P_voxels = V)
  expect_equal(res, expected, tolerance = 1e-12)

  perm_res <- compute_crossnobis_distances_sl(U_folds[, , sample(M)], P_voxels = V)
  expect_equal(perm_res, expected, tolerance = 1e-12)

  Iw <- diag(V)
  U_folds_I <- U_folds
  for (m in seq_len(M)) U_folds_I[, , m] <- U_folds_I[, , m] %*% Iw
  res_I <- compute_crossnobis_distances_sl(U_folds_I, P_voxels = V)
  expect_equal(res_I, expected, tolerance = 1e-12)
})


test_that("crossnobis distance is less biased than naive within-fold distance", {
  set.seed(456)
  K <- 2; V <- 5; M <- 4
  mu_true <- matrix(c(rep(0, V), rep(1, V)), nrow = 2, byrow = TRUE,
                    dimnames = list(paste0("C", 1:K), NULL))
  U_folds <- array(NA_real_, dim = c(K, V, M),
                   dimnames = list(paste0("C", 1:K), NULL, paste0("F", 1:M)))
  for (m in seq_len(M)) {
    U_folds[1, , m] <- mu_true[1, ] + rnorm(V, sd = 0.5)
    U_folds[2, , m] <- mu_true[2, ] + rnorm(V, sd = 0.5)
  }

  true_dist <- sum((mu_true[1, ] - mu_true[2, ])^2) / V
  delta <- U_folds[1, , ] - U_folds[2, , ]
  naive <- mean(colSums(delta^2)) / V
  cross_res <- compute_crossnobis_distances_sl(U_folds, P_voxels = V)["C1_vs_C2"]

  expect_lt(abs(cross_res - true_dist), abs(naive - true_dist))
})

