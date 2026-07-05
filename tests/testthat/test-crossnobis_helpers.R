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

# A more direct test for pair ordering matching lower.tri
test_that("compute_crossnobis_distances_sl returns pairs in lower.tri order", {
  K <- 3; V <- 1; M <- 2
  # create U_folds array filled with zeros
  U_folds <- array(0, dim = c(K, V, M),
                   dimnames = list(paste0("Cond", 1:K), NULL, NULL))
  res <- compute_crossnobis_distances_sl(U_folds, P_voxels = V)

  pair_idx <- which(lower.tri(matrix(1, K, K)), arr.ind = TRUE)
  expected_names <- paste0(paste0("Cond", pair_idx[,1]), "_vs_", paste0("Cond", pair_idx[,2]))
  expect_equal(names(res), expected_names)
  expect_equal(length(res), K * (K - 1) / 2)
})

# Helper to manually compute crossnobis distances
manual_crossnobis <- function(U_folds) {
  K <- dim(U_folds)[1]; V <- dim(U_folds)[2]; M <- dim(U_folds)[3]
  # Use lower.tri indices instead of combn to match the implementation
  pair_idx <- which(lower.tri(matrix(1, K, K)), arr.ind = TRUE)
  n_pairs <- nrow(pair_idx)
  res <- numeric(n_pairs)
  pair_names <- character(n_pairs)
  
  for (p in seq_len(n_pairs)) {
    i <- pair_idx[p, 1]; j <- pair_idx[p, 2]
    deltas <- U_folds[i, , ] - U_folds[j, , ]
    ip_mat <- crossprod(deltas)
    off_diag <- sum(ip_mat) - sum(diag(ip_mat))
    res[p] <- off_diag / (V * M * (M - 1))
    pair_names[p] <- paste0(dimnames(U_folds)[[1]][i], "_vs_", dimnames(U_folds)[[1]][j])
  }
  names(res) <- pair_names
  res
}

rsatoolbox_leave_one_out_reference <- function(U_folds, precision = NULL) {
  K <- dim(U_folds)[1]
  V <- dim(U_folds)[2]
  M <- dim(U_folds)[3]
  if (is.null(precision)) {
    precision <- diag(V)
  }

  pair_idx <- which(lower.tri(matrix(TRUE, nrow = K, ncol = K)), arr.ind = TRUE)
  rdms <- matrix(NA_real_, nrow = M, ncol = nrow(pair_idx))
  G_sum <- matrix(0, nrow = K, ncol = K,
                  dimnames = list(dimnames(U_folds)[[1]], dimnames(U_folds)[[1]]))

  for (m in seq_len(M)) {
    train <- apply(U_folds[, , -m, drop = FALSE], c(1, 2), mean)
    test <- U_folds[, , m, drop = FALSE][, , 1]
    kernel <- train %*% precision %*% t(test)
    dmat <- outer(diag(kernel), diag(kernel), "+") - kernel - t(kernel)
    rdms[m, ] <- dmat[lower.tri(dmat)] / V
    G_sum <- G_sum + kernel
  }

  pair_names <- paste0(dimnames(U_folds)[[1]][pair_idx[, 1]],
                       "_vs_",
                       dimnames(U_folds)[[1]][pair_idx[, 2]])
  distances <- colMeans(rdms)
  names(distances) <- pair_names

  list(
    distances = distances,
    second_moment = G_sum / M
  )
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

test_that("compute_crossnobis_second_moment_sl matches manual Gram and distances derive from it", {
  set.seed(321)
  K <- 4; V <- 3; M <- 4
  U_folds <- array(rnorm(K * V * M), dim = c(K, V, M),
                   dimnames = list(paste0("C", 1:K), paste0("V", 1:V), paste0("F", 1:M)))

  U_sum <- rowSums(U_folds, dims = 2)
  dim(U_sum) <- c(K, V)
  sum_gram <- tcrossprod(U_sum)
  within_fold_gram <- matrix(0, K, K)
  for (m in seq_len(M)) {
    within_fold_gram <- within_fold_gram + tcrossprod(U_folds[, , m])
  }
  expected_G <- (sum_gram - within_fold_gram) / (M * (M - 1))
  dimnames(expected_G) <- list(paste0("C", 1:K), paste0("C", 1:K))

  G_cv <- rMVPA:::compute_crossnobis_second_moment_sl(U_folds, vectorize = FALSE)
  expect_equal(G_cv, expected_G, tolerance = 1e-12)

  d_from_G <- (diag(G_cv)[row(G_cv)] + diag(G_cv)[col(G_cv)] - 2 * G_cv)[lower.tri(G_cv)] / V
  d_direct <- compute_crossnobis_distances_sl(U_folds, P_voxels = V)
  names(d_from_G) <- names(d_direct)
  expect_equal(d_direct, d_from_G, tolerance = 1e-12)
})

test_that("crossnobis helpers match rsatoolbox leave-one-out train/test formulation", {
  set.seed(224)
  K <- 5; V <- 4; M <- 3
  U_folds <- array(rnorm(K * V * M), dim = c(K, V, M),
                   dimnames = list(paste0("C", 1:K), paste0("V", 1:V), paste0("F", 1:M)))
  A <- matrix(rnorm(V * V), nrow = V)
  precision <- crossprod(A) + diag(V) * 0.25

  ref_euclidean <- rsatoolbox_leave_one_out_reference(U_folds)
  expect_equal(
    compute_crossnobis_distances_sl(U_folds, P_voxels = V),
    ref_euclidean$distances,
    tolerance = 1e-12
  )
  expect_equal(
    rMVPA:::compute_crossnobis_second_moment_sl(U_folds, vectorize = FALSE),
    ref_euclidean$second_moment,
    tolerance = 1e-12
  )

  W <- t(chol(precision))
  U_whitened <- U_folds
  for (m in seq_len(M)) {
    U_whitened[, , m] <- U_whitened[, , m] %*% W
  }

  ref_mahalanobis <- rsatoolbox_leave_one_out_reference(U_folds, precision = precision)
  expect_equal(
    compute_crossnobis_distances_sl(U_whitened, P_voxels = V),
    ref_mahalanobis$distances,
    tolerance = 1e-12
  )
  expect_equal(
    rMVPA:::compute_crossnobis_second_moment_sl(U_whitened, vectorize = FALSE),
    ref_mahalanobis$second_moment,
    tolerance = 1e-12
  )
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
  
  # Get the correct pair name based on lower.tri order
  pair_idx <- which(lower.tri(matrix(1, K, K)), arr.ind = TRUE)
  pair_name <- paste0("C", pair_idx[1, 1], "_vs_", "C", pair_idx[1, 2])
  cross_res <- compute_crossnobis_distances_sl(U_folds, P_voxels = V)[pair_name]

  expect_lt(abs(cross_res - true_dist), abs(naive - true_dist))
})
