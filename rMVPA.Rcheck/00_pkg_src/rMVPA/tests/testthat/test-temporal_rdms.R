library(testthat)
library(rMVPA)

context("Temporal RDMs for RSA nuisance modeling")

test_that("temporal_rdm creates valid dist objects", {
  # Test basic dist creation
  n <- 20
  idx <- 1:n
  
  temp_rdm <- temporal_rdm(idx, kernel = "adjacent", width = 1)
  
  expect_s3_class(temp_rdm, "dist")
  expect_equal(attr(temp_rdm, "Size"), n)
  expect_equal(length(temp_rdm), n * (n - 1) / 2)
})

test_that("temporal_rdm kernels produce expected patterns", {
  n <- 10
  idx <- 1:n
  
  # Test adjacent kernel
  temp_adj <- temporal_rdm(idx, kernel = "adjacent", width = 1, 
                           normalize = "none", as_dist = FALSE)
  expect_equal(temp_adj[1, 2], 1)  # Adjacent elements
  expect_equal(temp_adj[1, 3], 0)  # Non-adjacent elements
  expect_equal(diag(temp_adj), rep(0, n))  # Diagonal is zero
  
  # Test boxcar kernel
  temp_box <- temporal_rdm(idx, kernel = "boxcar", width = 2,
                           normalize = "none", as_dist = FALSE)
  expect_equal(temp_box[1, 2], 1)  # Within width
  expect_equal(temp_box[1, 3], 1)  # Within width
  expect_equal(temp_box[1, 4], 0)  # Outside width
  
  # Test linear kernel
  temp_lin <- temporal_rdm(idx, kernel = "linear",
                           normalize = "none", as_dist = FALSE)
  expect_equal(temp_lin[1, 2], 1)  # Lag of 1
  expect_equal(temp_lin[1, 5], 4)  # Lag of 4
  
  # Test exponential kernel
  temp_exp <- temporal_rdm(idx, kernel = "exp", lambda = 1,
                           normalize = "none", as_dist = FALSE)
  expect_equal(temp_exp[1, 2], exp(-1))  # Lag of 1
  expect_equal(temp_exp[1, 3], exp(-2))  # Lag of 2
})

test_that("temporal_rdm handles blocks correctly", {
  n <- 12
  idx <- 1:n
  blocks <- rep(1:3, each = 4)
  
  # Test within_blocks_only
  temp_within <- temporal_rdm(idx, block = blocks, kernel = "adjacent",
                              within_blocks_only = TRUE, 
                              normalize = "none", as_dist = FALSE)
  
  # Check that across-block pairs are zero
  expect_equal(temp_within[1, 5], 0)  # Block 1 to Block 2
  expect_equal(temp_within[4, 5], 0)  # Block 1 to Block 2
  
  # Check that within-block pairs are non-zero for adjacent
  expect_equal(temp_within[1, 2], 1)  # Within Block 1
  expect_equal(temp_within[5, 6], 1)  # Within Block 2
})

test_that("temporal_rdm normalization works", {
  n <- 10
  idx <- 1:n
  
  # Test rank normalization
  temp_rank <- temporal_rdm(idx, kernel = "exp", normalize = "rank", lambda = 0.5)
  vals_rank <- as.vector(temp_rank)
  # With exponential kernel, we should get unique values (no ties)
  
  # Test z-score normalization
  temp_z <- temporal_rdm(idx, kernel = "linear", normalize = "z")
  vals_z <- as.vector(temp_z)
  expect_equal(mean(vals_z), 0, tolerance = 1e-10)
  expect_equal(sd(vals_z), 1, tolerance = 1e-10)
})

test_that("temporal_rdm circular wrapping works", {
  n <- 8
  idx <- 1:n
  
  temp_wrap <- temporal_rdm(idx, kernel = "linear", wrap = TRUE,
                            normalize = "none", as_dist = FALSE)
  
  # With wrapping, distance from 1 to 8 should be 1 (circular)
  expect_equal(temp_wrap[1, 8], 1)
  # Distance from 1 to 5 should be min(4, 8-4) = 4
  expect_equal(temp_wrap[1, 5], 4)
})

test_that("temporal_nuisance_for_msreve produces correct dimensions", {
  # Create a simple mvpa_design mock
  n_samples <- 20
  n_conditions <- 4
  
  # Mock mvpa_design structure
  condition_labels <- factor(rep(paste0("cond", 1:n_conditions), 
                                 each = n_samples / n_conditions))
  blocks <- rep(1:2, each = n_samples / 2)
  
  mvpa_des <- list(
    Y = condition_labels,
    block_var = blocks,
    train_design = data.frame(cond = condition_labels, block = blocks),
    ncond = n_conditions
  )
  class(mvpa_des) <- "mvpa_design"
  
  time_idx <- 1:n_samples
  
  # Test basic functionality
  temp_K <- temporal_nuisance_for_msreve(mvpa_des, time_idx,
                                         reduce = "min", kernel = "exp")
  
  expect_true(is.matrix(temp_K))
  expect_equal(dim(temp_K), c(n_conditions, n_conditions))
  expect_equal(diag(temp_K), rep(0, n_conditions))
  expect_true(isSymmetric(temp_K))
})

test_that("temporal_nuisance_for_msreve reduction methods work", {
  n_samples <- 16
  n_conditions <- 4
  
  # Create structured conditions and time indices
  condition_labels <- factor(rep(paste0("cond", 1:n_conditions), 
                                 each = n_samples / n_conditions))
  blocks <- rep(1:2, each = n_samples / 2)
  time_idx <- 1:n_samples
  
  mvpa_des <- list(
    Y = condition_labels,
    block_var = blocks,
    train_design = data.frame(cond = condition_labels, block = blocks),
    ncond = n_conditions
  )
  class(mvpa_des) <- "mvpa_design"
  
  # Test different reduction methods
  temp_K_min <- temporal_nuisance_for_msreve(mvpa_des, time_idx,
                                             reduce = "min", kernel = "linear",
                                             within_blocks_only = FALSE)
  
  temp_K_mean <- temporal_nuisance_for_msreve(mvpa_des, time_idx,
                                              reduce = "mean", kernel = "linear",
                                              within_blocks_only = FALSE)
  
  # The mean should generally be >= min for temporal distances
  expect_true(all(temp_K_mean >= temp_K_min - 1e-10))
})

test_that("temporal wrapper function works", {
  n <- 10
  idx <- 1:n
  
  # Test that temporal() produces same result as temporal_rdm()
  temp1 <- temporal(idx, kernel = "exp", lambda = 2)
  temp2 <- temporal_rdm(idx, kernel = "exp", lambda = 2)
  
  expect_equal(as.vector(temp1), as.vector(temp2))
})

test_that("Integration with rsa_design works", {
  n <- 50
  blocks <- rep(1:5, each = 10)
  
  # Create a dummy task RDM
  task_rdm <- dist(matrix(rnorm(n * 10), n, 10))
  
  # Create temporal RDM
  temp_rdm <- temporal_rdm(1:n, block = blocks, kernel = "adjacent",
                           within_blocks_only = TRUE)
  
  # Create RSA design with temporal nuisance
  rdes <- rsa_design(~ task_rdm + temp_rdm,
                    data = list(task_rdm = task_rdm, 
                               temp_rdm = temp_rdm,
                               block = blocks),
                    block_var = "block",
                    keep_intra_run = TRUE)
  
  expect_s3_class(rdes, "rsa_design")
  expect_equal(length(rdes$model_mat), 2)
  expect_true("task_rdm" %in% names(rdes$model_mat))
  expect_true("temp_rdm" %in% names(rdes$model_mat))
})

test_that("Integration with msreve_design works", {
  # Create mock mvpa_design
  n_samples <- 16
  n_conditions <- 4
  
  condition_labels <- factor(rep(paste0("cond", 1:n_conditions), 
                                 each = n_samples / n_conditions))
  blocks <- rep(1:2, each = n_samples / 2)
  
  mvpa_des <- list(
    Y = condition_labels,
    block_var = blocks,
    train_design = data.frame(cond = condition_labels, block = blocks),
    ncond = n_conditions,
    conditions = levels(condition_labels)
  )
  class(mvpa_des) <- "mvpa_design"
  
  # Create contrast matrix
  C_mat <- matrix(c(1, 1, -1, -1,
                   1, -1, 0, 0), 
                 nrow = n_conditions, ncol = 2)
  rownames(C_mat) <- levels(condition_labels)
  colnames(C_mat) <- c("C1", "C2")
  
  # Create temporal nuisance
  temp_K <- temporal_nuisance_for_msreve(mvpa_des, 
                                         time_idx = 1:n_samples,
                                         reduce = "min",
                                         kernel = "exp", lambda = 2)
  
  # Create msreve_design with nuisance
  msreve_des <- msreve_design(mvpa_des, C_mat,
                              nuisance_rdms = list(temporal = temp_K))
  
  expect_s3_class(msreve_des, "msreve_design")
  expect_true(!is.null(msreve_des$nuisance_rdms))
  expect_equal(names(msreve_des$nuisance_rdms), "temporal")
  expect_equal(dim(msreve_des$nuisance_rdms$temporal), c(n_conditions, n_conditions))
})

test_that("Edge cases are handled correctly", {
  # Test with single observation
  expect_error(temporal_rdm(1), NA)  # Should not error
  
  # Test with two observations
  temp2 <- temporal_rdm(1:2, kernel = "linear", as_dist = TRUE)
  expect_equal(length(temp2), 1)
  
  # Test with all same block
  n <- 10
  idx <- 1:n
  blocks <- rep(1, n)
  temp_same <- temporal_rdm(idx, block = blocks, within_blocks_only = TRUE,
                            kernel = "adjacent")
  expect_true(all(as.vector(temp_same) >= 0))  # All within same block
  
  # Test with all different blocks (each sample in its own block)
  blocks_diff <- 1:n
  temp_diff <- temporal_rdm(idx, block = blocks_diff, within_blocks_only = TRUE,
                           kernel = "adjacent", normalize = "none", as_dist = FALSE)
  # When each sample is in its own block and within_blocks_only = TRUE,
  # only diagonal elements are within-block (which are always 0)
  # So all off-diagonal elements should be 0
  diag(temp_diff) <- NA  # Ignore diagonal
  expect_true(all(temp_diff == 0, na.rm = TRUE))  # All off-diagonal are 0
})

test_that("Parameter validation works", {
  n <- 10
  idx <- 1:n
  
  # Test invalid kernel
  expect_error(temporal_rdm(idx, kernel = "invalid"))
  
  # Test invalid normalization
  expect_error(temporal_rdm(idx, normalize = "invalid"))
  
  # Test mismatched block length
  expect_error(temporal_rdm(idx, block = rep(1, n-1)))
})

test_that("Nuisance RDMs in contrast_rsa_model don't affect output dimensions", {
  skip_if_not_installed("neuroim2")
  
  # Create minimal test data
  set.seed(123)
  n_samples <- 16
  n_voxels <- 10
  n_conditions <- 4
  
  # Create dummy data as NeuroVec so mvpa_dataset accepts it
  dummy_array <- array(rnorm(n_voxels * n_samples), c(n_voxels, 1, 1, n_samples))
  dummy_space <- neuroim2::NeuroSpace(c(n_voxels, 1, 1, n_samples))
  dummy_sl_vec <- neuroim2::NeuroVec(dummy_array, dummy_space)
  dummy_mask <- neuroim2::NeuroVol(array(1, c(n_voxels, 1, 1)), neuroim2::NeuroSpace(c(n_voxels, 1, 1)))
  
  condition_labels <- factor(rep(paste0("cond", 1:n_conditions), 
                                 each = n_samples / n_conditions))
  run_labels <- factor(rep(1:2, each = n_samples / 2))
  
  # Create mvpa objects
  mvpa_dat <- mvpa_dataset(train_data = dummy_sl_vec, mask = dummy_mask)
  mvpa_des <- mvpa_design(
    train_design = data.frame(condition = condition_labels, run = run_labels),
    y_train = ~ condition,
    block_var = ~ run
  )
  
  # Create contrast matrix
  C_mat <- matrix(c(1, 1, -1, -1), nrow = n_conditions, ncol = 1)
  condition_levels <- levels(rMVPA::y_train(mvpa_des))
  rownames(C_mat) <- condition_levels
  colnames(C_mat) <- "Contrast1"
  
  # Test without nuisance
  msreve_des1 <- msreve_design(mvpa_des, C_mat)
  model1 <- contrast_rsa_model(mvpa_dat, msreve_des1, 
                               output_metric = "beta_only")
  
  # Test with nuisance
  temp_K <- matrix(rnorm(n_conditions^2), n_conditions, n_conditions)
  temp_K <- (temp_K + t(temp_K)) / 2  # Make symmetric
  diag(temp_K) <- 0
  
  msreve_des2 <- msreve_design(mvpa_des, C_mat,
                               nuisance_rdms = list(temp = temp_K))
  model2 <- contrast_rsa_model(mvpa_dat, msreve_des2,
                               output_metric = "beta_only")
  
  # Both models should have same structure
  expect_equal(class(model1), class(model2))
  expect_equal(length(model1$output_metric), length(model2$output_metric))
})
