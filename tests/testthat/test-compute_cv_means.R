context("compute_crossvalidated_means_sl")

library(testthat)
library(rMVPA)
library(tibble)
library(assertthat)

# --- Helper objects for tests ---

# Mock mvpa_design (simplified for this test)
mock_mvpa_design_cv <- function(n_samples = 24, n_cond = 4, n_blocks = 3, Y_labels=NULL) {
  if (is.null(Y_labels)) {
    cond_labels <- factor(rep(paste0("Cond", 1:n_cond), each = n_samples / n_cond))
  } else {
    cond_labels <- factor(Y_labels)
    n_cond <- length(levels(cond_labels))
  }
  block_var <- factor(rep(1:n_blocks, length.out = n_samples)) # ensure block_var has n_samples length
  conditions <- levels(cond_labels)
  # design_matrix is used for nrow check in the function
  design_matrix <- matrix(0, nrow = n_samples, ncol = 1) 
  
  structure(
    list(
      Y = cond_labels,
      block_var = block_var,
      conditions = conditions,
      ncond = n_cond,
      design_matrix = design_matrix, # Added for input check
      samples = 1:n_samples
    ),
    class = c("mvpa_design", "list")
  )
}

# Mock cv_spec object with S3 methods
mock_cv_spec_s3 <- function(mvpa_design) {
  n_folds_val <- length(unique(mvpa_design$block_var))
  folds_val <- as.integer(mvpa_design$block_var)
  
  obj <- list(
    # Store necessary data for methods
    .folds_val = folds_val,
    .n_folds_val = n_folds_val
  )
  class(obj) <- c("mock_cv_spec", "cross_validation", "list")
  obj
}

#' @export
get_nfolds.mock_cv_spec <- function(obj, ...) {
  obj$.n_folds_val
}

#' @export
train_indices.mock_cv_spec <- function(obj, fold_num, ...) {
  which(obj$.folds_val != fold_num)
}

# --- Tests for compute_crossvalidated_means_sl ---

test_that("compute_crossvalidated_means_sl works with 'average' method and return_folds=FALSE", {
  n_samples <- 12
  n_cond <- 3
  n_blocks <- 3 # 3 blocks, 4 samples per block
  n_voxels <- 5
  
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  
  set.seed(123)
  # Each condition has a base mean (1, 2, or 3) across all voxels
  # Block 1: Cond1, Cond1, Cond1, Cond1 (all get base 1)
  # Block 2: Cond2, Cond2, Cond2, Cond2 (all get base 2)
  # Block 3: Cond3, Cond3, Cond3, Cond3 (all get base 3)
  # This setup is not ideal as a condition is only in one block.
  # Let's adjust sample and block structure for better CV test
  n_samples <- 12; n_cond <- 3; n_blocks <- 2 # 2 blocks, 6 obs/block, 2 obs/cond/block
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)

  cond_means_per_sample <- as.numeric(mvpa_des$Y) # Cond1=1, Cond2=2, Cond3=3
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), nrow = n_samples, ncol = n_voxels) + 
               matrix(rnorm(n_samples * n_voxels, 0, 0.01), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)
  
  U_hat <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = FALSE)
  
  expect_equal(dim(U_hat), c(n_cond, n_voxels))
  expect_equal(rownames(U_hat), paste0("Cond", 1:n_cond))
  # Check that U_hat for CondX is close to X (e.g. U_hat[1,] should be ~1)
  expect_true(all(abs(U_hat[1,] - 1) < 0.1))
  expect_true(all(abs(U_hat[2,] - 2) < 0.1))
  expect_true(all(abs(U_hat[3,] - 3) < 0.1))
})

test_that("compute_crossvalidated_means_sl works with 'average' method and return_folds=TRUE", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)

  set.seed(1234)
  cond_means_per_sample <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), nrow = n_samples, ncol = n_voxels) + 
               matrix(rnorm(n_samples * n_voxels, 0, 0.01), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)

  res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
  
  expect_true(is.list(res_list))
  expect_equal(names(res_list), c("mean_estimate", "fold_estimates"))
  
  U_hat_mean <- res_list$mean_estimate
  U_folds <- res_list$fold_estimates
  
  expect_equal(dim(U_hat_mean), c(n_cond, n_voxels))
  expect_equal(rownames(U_hat_mean), paste0("Cond", 1:n_cond))
  expect_true(all(abs(U_hat_mean[1,] - 1) < 0.1))
  expect_true(all(abs(U_hat_mean[2,] - 2) < 0.1))
  expect_true(all(abs(U_hat_mean[3,] - 3) < 0.1))
  
  expect_equal(dim(U_folds), c(n_cond, n_voxels, n_blocks))
  expect_equal(dimnames(U_folds)[[1]], paste0("Cond", 1:n_cond))
  expect_equal(dimnames(U_folds)[[2]], paste0("V", 1:n_voxels))
  expect_equal(dimnames(U_folds)[[3]], paste0("Fold", 1:n_blocks))
  
  # Check if mean_estimate is roughly the mean of fold_estimates, handling NAs
  # This is complex due to the online averaging. A simpler check:
  # For each fold, the training data for CondX should average to X.
  # For fold 1 (train on block 2): block_var for block2 is 2
  # Samples in block 2: mvpa_des$Y[mvpa_des$block_var==2]
  # Cond1 samples in block 2: sl_data[mvpa_des$Y=="Cond1" & mvpa_des$block_var==2, ]
  # This value should be reflected in U_folds["Cond1", , 1]
  # This is also complex. For now, basic structure checks are primary.
  expect_false(any(is.na(U_folds))) # With this balanced design, no NAs expected
})

test_that("compute_crossvalidated_means_sl with return_folds=TRUE matches return_folds=FALSE for mean_estimate", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(125)
  cond_means_per_sample <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), nrow = n_samples, ncol = n_voxels) + 
               matrix(rnorm(n_samples * n_voxels, 0, 0.01), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)
  
  U_hat_direct <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = FALSE)
  res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
  
  expect_equal(U_hat_direct, res_list$mean_estimate)
})


test_that("compute_crossvalidated_means_sl handles NAs correctly with return_folds=TRUE", {
  n_samples <- 6 
  n_cond <- 3
  n_blocks <- 2
  n_voxels <- 2
  
  # Cond1 in B1, Cond2 in B1, Cond3 in B2
  # Training on B2 means C1, C2 are NA. Training on B1 means C3 is NA.
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks, Y_labels=c("C1", "C2", "C1", "C2", "C3", "C3"))
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  
  set.seed(126)
  sl_data <- matrix(rnorm(n_samples * n_voxels), nrow = n_samples, ncol=n_voxels)
  colnames(sl_data) <- paste0("V",1:n_voxels)

  res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
  
  U_hat_mean <- res_list$mean_estimate
  U_folds <- res_list$fold_estimates

  # Given Y_labels=c("C1", "C2", "C1", "C2", "C3", "C3") and block_var=c(1,2,1,2,1,2)
  # Fold 1 trains on block 2 data (samples 2,4,6 with labels C2,C2,C3)
  # Fold 2 trains on block 1 data (samples 1,3,5 with labels C1,C1,C3)

  # U_hat_mean expectations:
  # C1: present in Fold 2 train, missing Fold 1 train. count=1. Should be NA.
  # C2: present in Fold 1 train, missing Fold 2 train. count=1. Should be NA.
  # C3: present in Fold 1 train AND Fold 2 train. count=2. Should NOT be NA.
  expect_true(all(is.na(U_hat_mean["C1",])))
  expect_true(all(is.na(U_hat_mean["C2",])))
  expect_false(all(is.na(U_hat_mean["C3",]))) # C3 is in both training folds

  # U_folds expectations for Fold 1 (trains on C2,C2,C3):
  expect_true(all(is.na(U_folds["C1", , 1])))    # C1 not in train
  expect_false(all(is.na(U_folds["C2", , 1])))   # C2 is in train
  expect_false(any(is.na(U_folds["C2", , 1])))    # C2 should be fully populated
  expect_false(all(is.na(U_folds["C3", , 1])))   # C3 is in train
  expect_false(any(is.na(U_folds["C3", , 1])))    # C3 should be fully populated
  
  # U_folds expectations for Fold 2 (trains on C1,C1,C3):
  expect_false(all(is.na(U_folds["C1", , 2])))   # C1 is in train
  expect_false(any(is.na(U_folds["C1", , 2])))    # C1 should be fully populated
  expect_true(all(is.na(U_folds["C2", , 2])))    # C2 not in train
  expect_false(all(is.na(U_folds["C3", , 2])))   # C3 is in train
  expect_false(any(is.na(U_folds["C3", , 2])))    # C3 should be fully populated
})

test_that("compute_crossvalidated_means_sl with 'crossnobis' and W, return_folds=TRUE", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 4 # V must be even for simple W
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(127)
  cond_means_per_sample <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), nrow = n_samples, ncol = n_voxels) + 
               matrix(rnorm(n_samples * n_voxels, 0, 0.05), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)

  # Simple whitening matrix (e.g. swap pairs of voxels)
  W <- diag(n_voxels)
  W[1,1] <- 0; W[2,2] <- 0; W[1,2] <- 1; W[2,1] <- 1 # Swap V1 and V2
  W[3,3] <- 0; W[4,4] <- 0; W[3,4] <- 1; W[4,3] <- 1 # Swap V3 and V4
  
  # Case 1: Crossnobis with W
  res_list_W <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                estimation_method = "crossnobis", 
                                                whitening_matrix_W = W, 
                                                return_folds = TRUE)
  expect_true(is.list(res_list_W))
  expect_equal(names(res_list_W), c("mean_estimate", "fold_estimates"))
  expect_equal(dim(res_list_W$mean_estimate), c(n_cond, n_voxels))
  expect_equal(dim(res_list_W$fold_estimates), c(n_cond, n_voxels, n_blocks))
  
  # Check that whitening had an effect (compare to average without W on fold_estimates)
  res_list_noW_folds_only <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                estimation_method = "average", 
                                                return_folds = TRUE)$fold_estimates
  
  # Whitening is applied per-fold. The res_list_W$fold_estimates should be different.
  expect_false(isTRUE(all.equal(res_list_W$fold_estimates, res_list_noW_folds_only)))
  
  # Example: U_folds_noW_f1_c1 = res_list_noW_folds_only["Cond1", , 1]
  # U_folds_W_f1_c1_expected = U_folds_noW_f1_c1 %*% W
  # This should hold for res_list_W$fold_estimates["Cond1", , 1]
  cond_idx <- 1; fold_idx <- 1
  unwhitened_fold_mean_cond <- res_list_noW_folds_only[cond_idx, , fold_idx]
  expected_whitened_fold_mean_cond <- matrix(unwhitened_fold_mean_cond, nrow=1) %*% W
  
  expect_equal(as.vector(res_list_W$fold_estimates[cond_idx, , fold_idx]), 
               as.vector(expected_whitened_fold_mean_cond), tolerance=1e-9)

  # Verify U_hat_mean is the mean of the (whitened) fold_estimates
  calculated_mean_from_folds <- apply(res_list_W$fold_estimates, c(1, 2), mean, na.rm = TRUE)
  # apply with mean, na.rm=TRUE can produce NaN if all inputs for a mean are NA. 
  # The function compute_crossvalidated_means_sl produces NA_real_ in such cases.
  calculated_mean_from_folds[is.nan(calculated_mean_from_folds)] <- NA_real_
  expect_equal(res_list_W$mean_estimate, calculated_mean_from_folds, tolerance=1e-1)

  # Case 2: Crossnobis expects W, should error if W is NULL (current behavior)
  expect_error(
    compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                    estimation_method = "crossnobis", 
                                    whitening_matrix_W = NULL, 
                                    return_folds = TRUE),
    "`whitening_matrix_W` must be provided when `estimation_method = crossnobis`"
  )
})

test_that("compute_crossvalidated_means_sl ('crossnobis') correctly whitens per-fold estimates voxel-wise", {
  n_samples <- 12
  n_cond <- 2 # Simpler with fewer conditions for this specific test
  n_blocks <- 2
  n_voxels <- 4 
  
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  
  set.seed(456) # New seed
  
  # Create sl_data with more distinct voxel patterns per condition
  sl_data <- matrix(0, nrow = n_samples, ncol = n_voxels)
  colnames(sl_data) <- paste0("V", 1:n_voxels)
  
  # Condition 1: Stronger on V1 & V3 initially
  # Condition 2: Stronger on V2 & V4 initially
  base_pattern_cond1 <- c(5, 1, 4, 2) 
  base_pattern_cond2 <- c(1, 6, 2, 5)
  
  for (i in 1:n_samples) {
    noise <- rnorm(n_voxels, 0, 0.1)
    if (mvpa_des$Y[i] == "Cond1") {
      sl_data[i, ] <- base_pattern_cond1 + noise
    } else { # Cond2
      sl_data[i, ] <- base_pattern_cond2 + noise
    }
  }

  # Whitening matrix W (swaps V1-V2 and V3-V4)
  W <- diag(n_voxels)
  W[1,1]<-0; W[2,2]<-0; W[1,2]<-1; W[2,1]<-1
  W[3,3]<-0; W[4,4]<-0; W[3,4]<-1; W[4,3]<-1
  
  # Get fold estimates from "crossnobis" (whitened)
  res_crossnobis <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                    estimation_method = "crossnobis", 
                                                    whitening_matrix_W = W, 
                                                    return_folds = TRUE)
  whitened_fold_estimates <- res_crossnobis$fold_estimates
  
  # Get fold estimates from "average" (unwhitened)
  res_average <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                 estimation_method = "average", 
                                                 return_folds = TRUE)
  unwhitened_fold_estimates <- res_average$fold_estimates
  
  expect_equal(dim(whitened_fold_estimates), c(n_cond, n_voxels, n_blocks))
  expect_equal(dim(unwhitened_fold_estimates), c(n_cond, n_voxels, n_blocks))

  # Iterate through conditions and folds to check whitening
  for (k in 1:n_cond) { # k for condition index
    cond_name <- rownames(unwhitened_fold_estimates)[k]
    for (f in 1:n_blocks) { # f for fold index
      
      # Get the unwhitened mean pattern for this condition (k) and fold (f)
      U_kf_unwhitened <- unwhitened_fold_estimates[cond_name, , f, drop = TRUE]
      
      # Manually whiten it
      # Ensure U_kf_unwhitened is a row vector for matrix multiplication
      U_kf_manual_whitened <- matrix(U_kf_unwhitened, nrow = 1) %*% W
      
      # Get the corresponding whitened pattern from the function's output
      U_kf_function_whitened <- whitened_fold_estimates[cond_name, , f, drop = TRUE]
      
      test_label <- paste0("Cond: ", cond_name, ", Fold: ", f)
      expect_equal(as.vector(U_kf_function_whitened), 
                   as.vector(U_kf_manual_whitened), 
                   tolerance = 1e-9, 
                   label = paste(test_label, "- whitened vs manual"))
    }
  }
})

# Test for 'L2_norm' method with return_folds=TRUE
test_that("compute_crossvalidated_means_sl works with 'L2_norm' method and return_folds=TRUE", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)

  set.seed(12345)
  cond_means_per_sample <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), nrow = n_samples, ncol = n_voxels) + 
               matrix(rnorm(n_samples * n_voxels, 0, 0.01), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)

  res_list_l2 <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "L2_norm", return_folds = TRUE)
  
  expect_true(is.list(res_list_l2))
  expect_equal(names(res_list_l2), c("mean_estimate", "fold_estimates"))
  
  U_hat_mean_l2 <- res_list_l2$mean_estimate
  U_folds_l2 <- res_list_l2$fold_estimates
  
  expect_equal(dim(U_hat_mean_l2), c(n_cond, n_voxels))
  # Check L2 norm of rows of mean_estimate
  row_norms_mean <- sqrt(rowSums(U_hat_mean_l2^2))
  expect_true(all(abs(row_norms_mean - 1) < 1e-9))
  
  # Fold estimates should NOT be L2 normalized by this specific call
  # (L2_norm applies to the final mean estimate)
  res_list_avg <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
  expect_equal(U_folds_l2, res_list_avg$fold_estimates)
  
  # Verify U_hat_mean_l2 is the L2 normalized version of the mean of U_folds_l2
  mean_of_raw_folds <- apply(U_folds_l2, c(1,2), mean, na.rm = TRUE)
  mean_of_raw_folds[is.nan(mean_of_raw_folds)] <- NA_real_ # Consistent NA handling
  
  # L2 normalize this mean (handle rows that were all NA before summing)
  is_row_all_na <- apply(mean_of_raw_folds, 1, function(row) all(is.na(row)))
  
  norm_denom <- sqrt(rowSums(mean_of_raw_folds^2, na.rm=TRUE))
  # Avoid division by zero or NA/NaN issues for rows that were all NA or zero sum
  norm_denom[norm_denom < .Machine$double.eps & !is_row_all_na] <- 1 
  
  expected_l2_mean_estimate <- mean_of_raw_folds / norm_denom
  # If a row in mean_of_raw_folds was all NA, division by norm_denom (which would be 0 from na.rm=T, then 1)
  # would make it numeric. Re-apply NAs for rows that were entirely NA.
  expected_l2_mean_estimate[is_row_all_na, ] <- NA_real_
  
  expect_equal(U_hat_mean_l2, expected_l2_mean_estimate, tolerance=1e-2)
})


# Original tests (adapted for new mock_cv_spec_s3 and explicit return_folds=FALSE)
test_that("Original: compute_crossvalidated_means_sl handles single voxel data", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 1
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(124)
  cond_means <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(cond_means, nrow = n_samples, ncol = n_voxels) + 
               matrix(rnorm(n_samples * n_voxels, 0, 0.1), nrow = n_samples)
  colnames(sl_data) <- "V1"
               
  U_hat <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds=FALSE)
  
  expect_equal(dim(U_hat), c(n_cond, n_voxels))
  expect_equal(rownames(U_hat), paste0("Cond", 1:n_cond))
  expect_true(all(abs(U_hat[,1] - c(1,2,3)) < 0.2)) # Wider tolerance for single voxel
})

test_that("Original: compute_crossvalidated_means_sl returns NAs if condition missing in training fold", {
  n_samples <- 6 
  n_cond <- 3
  n_blocks <- 2
  n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks, Y_labels=c("C1", "C2", "C1", "C2", "C3", "C3"))
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(125)
  sl_data <- matrix(rnorm(n_samples * n_voxels), nrow = n_samples)
  colnames(sl_data) <- paste0("V",1:n_voxels)
  
  U_hat <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds=FALSE)

  # Expectations based on Y_labels=c("C1", "C2", "C1", "C2", "C3", "C3") and block_var=c(1,2,1,2,1,2)
  # C1: present in Fold 2 train, missing Fold 1 train. count=1. Final U_hat should be NA.
  # C2: present in Fold 1 train, missing Fold 2 train. count=1. Final U_hat should be NA.
  # C3: present in Fold 1 train AND Fold 2 train. count=2. Final U_hat should NOT be NA.
  expect_equal(dim(U_hat), c(n_cond, n_voxels))
  expect_true(all(is.na(U_hat[1,])))
  expect_true(all(is.na(U_hat[2,])))
  expect_false(all(is.na(U_hat[3,]))) # Corrected: C3 is in both training folds
})


test_that("Original: compute_crossvalidated_means_sl errors on invalid inputs", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  sl_data <- matrix(rnorm(n_samples * n_voxels), nrow = n_samples)
  
  expect_error(compute_crossvalidated_means_sl(as.data.frame(sl_data), mvpa_des, cv_spec, "average",return_folds=FALSE))
  # design_matrix check is now in mvpa_design mock
  mvpa_des_bad_rows <- mvpa_des
  mvpa_des_bad_rows$design_matrix <- matrix(0, nrow=n_samples-1, ncol=1)
  expect_error(compute_crossvalidated_means_sl(sl_data, mvpa_des_bad_rows, cv_spec, "average", return_folds=FALSE))
  
  expect_error(compute_crossvalidated_means_sl(sl_data, list(), cv_spec, "average", return_folds=FALSE))
  expect_error(compute_crossvalidated_means_sl(sl_data, mvpa_des, list(), "average", return_folds=FALSE))
  expect_error(compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "wrong_method", return_folds=FALSE))
  expect_error(compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = "not_logical"))
}) 