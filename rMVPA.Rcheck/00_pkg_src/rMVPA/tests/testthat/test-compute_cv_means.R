context("compute_crossvalidated_means_sl")

library(testthat)
library(rMVPA)
library(tibble)
library(assertthat)

# --- Helper objects for tests ---

# Mock mvpa_design (simplified for this test)
mock_mvpa_design_cv <- function(n_samples = 24, n_cond = 4, n_blocks = 3, Y_labels=NULL, custom_block_var=NULL) {
  if (is.null(Y_labels)) {
    cond_labels <- factor(rep(paste0("Cond", 1:n_cond), each = n_samples / n_cond))
  } else {
    cond_labels <- factor(Y_labels)
    n_cond <- length(levels(cond_labels))
  }
  if (is.null(custom_block_var)) {
    block_var <- factor(rep(1:n_blocks, length.out = n_samples))
  } else {
    block_var <- factor(custom_block_var)
    n_blocks <- length(unique(block_var)) # Update n_blocks if custom_block_var is provided
  }
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

# --- Tests for compute_crossvalidated_means_sl ---

# Mock implementations for with_mocked_bindings
.mock_get_nfolds <- function(obj, ...) {
  if (inherits(obj, "mock_cv_spec")) {
    return(obj$.n_folds_val)
  }
  stop(".mock_get_nfolds called with unexpected object type in this test context.")
}

.mock_train_indices <- function(obj, fold_num, ...) {
  if (inherits(obj, "mock_cv_spec")) {
    return(which(obj$.folds_val != fold_num))
  }
  stop(".mock_train_indices called with unexpected object type in this test context.")
}

test_that("compute_crossvalidated_means_sl works with 'average' method and return_folds=FALSE", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des) 
  
  set.seed(123)
  cond_means_per_sample <- as.numeric(mvpa_des$Y) 
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), nrow = n_samples, ncol = n_voxels) + 
               matrix(rnorm(n_samples * n_voxels, 0, 0.01), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)
  
  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      U_hat <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = FALSE)
      
      expect_equal(dim(U_hat), c(n_cond, n_voxels))
      expect_equal(rownames(U_hat), paste0("Cond", 1:n_cond))
      expect_true(all(abs(U_hat[1,] - 1) < 0.1))
      expect_true(all(abs(U_hat[2,] - 2) < 0.1))
      expect_true(all(abs(U_hat[3,] - 3) < 0.1))
    }
  )
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

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
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
      expect_false(any(is.na(U_folds)))
    }
  )
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

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      U_hat_direct <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = FALSE)
      res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
      expect_equal(U_hat_direct, res_list$mean_estimate)
    }
  )
})

test_that("process_single_fold reproduces fold_estimates from compute_crossvalidated_means_sl", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 4
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(999)
  cond_means <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(rep(cond_means, n_voxels), nrow = n_samples, ncol = n_voxels) +
              matrix(rnorm(n_samples * n_voxels, 0, 0.01), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      res <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
      for (f in seq_len(.mock_get_nfolds(cv_spec))) {
        train_idx <- .mock_train_indices(cv_spec, f)
        fold_expected <- rMVPA:::process_single_fold(sl_data[train_idx, , drop = FALSE],
                                                     mvpa_des$Y[train_idx],
                                                     levels(mvpa_des$Y),
                                                     length(levels(mvpa_des$Y)),
                                                     ncol(sl_data),
                                                     "average",
                                                     NULL,
                                                     colnames(sl_data))
        expect_equal(res$fold_estimates[,,f], fold_expected)
      }
    }
  )
})


test_that("compute_crossvalidated_means_sl handles NAs correctly with return_folds=TRUE", {
  n_samples <- 6 
  n_cond <- 3
  n_blocks <- 2
  n_voxels <- 2
  
  # Y_labels provide "C1", "C2", "C3" as condition names
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks, Y_labels=c("C1", "C2", "C1", "C2", "C3", "C3"))
  cv_spec <- mock_cv_spec_s3(mvpa_des) # block_var becomes c(1,2,1,2,1,2)
  
  set.seed(126)
  sl_data <- matrix(rnorm(n_samples * n_voxels), nrow = n_samples, ncol=n_voxels)
  colnames(sl_data) <- paste0("V",1:n_voxels)

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
      
      U_hat_mean <- res_list$mean_estimate
      U_folds <- res_list$fold_estimates

      # Fold 1 trains on samples with block_var == 2 (samples 2, 4, 6 with Y labels "C2", "C2", "C3")
      expect_true(all(is.na(U_folds["C1", , 1]))) # C1 not in training for Fold 1
      expect_false(any(is.na(U_folds["C2", , 1])))# C2 is in training for Fold 1
      expect_false(any(is.na(U_folds["C3", , 1])))# C3 is in training for Fold 1

      # Fold 2 trains on samples with block_var == 1 (samples 1, 3, 5 with Y labels "C1", "C1", "C3")
      expect_false(any(is.na(U_folds["C1", , 2])))# C1 is in training for Fold 2
      expect_true(all(is.na(U_folds["C2", , 2]))) # C2 not in training for Fold 2
      expect_false(any(is.na(U_folds["C3", , 2])))# C3 is in training for Fold 2

      # Mean estimate logic for 'average' method:
      # NA if cond_fold_counts < n_folds (which is 2)
      # C1: Present in Fold 2 training (count=1). Should be NA.
      # C2: Present in Fold 1 training (count=1). Should be NA.
      # C3: Present in Fold 1 & Fold 2 training (count=2). Should NOT be NA.
      expect_true(all(is.na(U_hat_mean["C1",])))
      expect_true(all(is.na(U_hat_mean["C2",])))
      expect_false(any(is.na(U_hat_mean["C3",])))
    }
  )
})

test_that("compute_crossvalidated_means_sl with 'crossnobis' and W, return_folds=TRUE", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(127)
  cond_means_per_sample <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), nrow = n_samples, ncol = n_voxels) +
               matrix(rnorm(n_samples * n_voxels, 0, 0.1), nrow = n_samples)
  colnames(sl_data) <- paste0("V", 1:n_voxels)
  W <- diag(n_voxels)

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                estimation_method = "crossnobis", 
                                                whitening_matrix_W = W, 
                                                return_folds = TRUE)
      
      expect_true(is.list(res_list))
      expect_equal(names(res_list), c("mean_estimate", "fold_estimates"))
      U_hat_mean <- res_list$mean_estimate
      U_folds <- res_list$fold_estimates
      
      expect_equal(dim(U_hat_mean), c(n_cond, n_voxels))
      expect_equal(dim(U_folds), c(n_cond, n_voxels, n_blocks))
      expect_false(any(is.na(U_folds)))
      expect_false(any(is.na(U_hat_mean)))

      # Test if whitening was applied (even if identity W in this specific setup)
      # For fold 1 (trains on block 2 data)
      train_idx_f1 <- .mock_train_indices(cv_spec, 1)
      # Get what the raw (unwhitened) means for fold 1 would have been
      raw_means_f1 <- rMVPA:::process_single_fold(sl_data[train_idx_f1, , drop=FALSE],
                                                  mvpa_des$Y[train_idx_f1],
                                                  levels(mvpa_des$Y), n_cond, n_voxels,
                                                  "average", NULL, colnames(sl_data))
      # Expected whitened means for fold 1 = raw_means_f1 %*% W
      expected_whitened_fold1 <- raw_means_f1 %*% W
      # Ensure colnames match for expect_equal attribute check
      colnames(expected_whitened_fold1) <- colnames(raw_means_f1) 
      
      expect_equal(U_folds[,,1], expected_whitened_fold1, tolerance = 1e-6)
    }
  )
})

test_that("compute_crossvalidated_means_sl correctly handles missing conditions in folds for 'average'", {
  n_samples <- 4; n_cond <- 2; n_blocks <- 2; n_voxels <- 2
  # Y_labels: C1, C1, C2, C2
  # block_var:1,  2,  1,  2 (samples 1,2,3,4)
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks, Y_labels=c("C1","C1","C2","C2"))
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  sl_data <- matrix(1:(n_samples*n_voxels), n_samples, n_voxels)

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
      U_hat_mean <- res_list$mean_estimate
      U_folds <- res_list$fold_estimates
      
      # Fold 1 trains on block 2 data (samples 2 and 4: Y labels C1, C2)
      expect_false(any(is.na(U_folds["C1", , 1]))) # C1 is in training for fold 1
      expect_false(any(is.na(U_folds["C2", , 1]))) # C2 is in training for fold 1
      
      # Fold 2 trains on block 1 data (samples 1 and 3: Y labels C1, C2)
      expect_false(any(is.na(U_folds["C1", , 2]))) # C1 is in training for fold 2
      expect_false(any(is.na(U_folds["C2", , 2]))) # C2 is in training for fold 2

      # Final U_hat_mean for 'average' method: NA if cond_fold_counts < n_folds.
      # Both C1 and C2 are in training for both folds (n_folds = 2).
      # So, cond_fold_counts for C1 and C2 should be 2. Neither should be NA.
      expect_false(any(is.na(U_hat_mean["C1",]))) 
      expect_false(any(is.na(U_hat_mean["C2",]))) 
    }
  )
})

test_that("compute_crossvalidated_means_sl ('crossnobis') correctly whitens per-fold estimates voxel-wise", {
  n_samples <- 12; n_cond <- 2; n_blocks <- 3; n_voxels <- 4
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(128)
  sl_data <- matrix(rnorm(n_samples * n_voxels), n_samples, n_voxels)
  W <- matrix(c(0,1,0,0, 1,0,0,0, 0,0,2,0, 0,0,0,0.5), n_voxels, n_voxels, byrow=TRUE)
  
  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                estimation_method = "crossnobis", 
                                                whitening_matrix_W = W, 
                                                return_folds = TRUE)
      U_folds_whitened <- res_list$fold_estimates
      
      for (f in 1:.mock_get_nfolds(cv_spec)) {
        train_idx <- .mock_train_indices(cv_spec, f)
        raw_fold_means <- rMVPA:::process_single_fold(sl_data[train_idx, , drop = FALSE],
                                                      mvpa_des$Y[train_idx],
                                                      levels(mvpa_des$Y), n_cond, n_voxels,
                                                      "average", NULL, colnames(sl_data))
        expected_whitened_fold_means <- raw_fold_means %*% W
        expect_equal(U_folds_whitened[,,f], expected_whitened_fold_means, tolerance = 1e-9)
      }
    }
  )
})

test_that("compute_crossvalidated_means_sl works with 'L2_norm' method and return_folds=TRUE", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 2; n_voxels <- 5
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  set.seed(129)
  cond_means_per_sample <- as.numeric(mvpa_des$Y)
  sl_data <- matrix(rep(cond_means_per_sample, n_voxels), n_samples, n_voxels) +
             matrix(rnorm(n_samples * n_voxels, 0, 0.1), n_samples, n_voxels)
  sl_data <- sl_data * sample(c(1,5,10), n_samples, replace=TRUE)

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                estimation_method = "L2_norm", 
                                                return_folds = TRUE)
      U_hat_norm <- res_list$mean_estimate
      U_folds_raw <- res_list$fold_estimates

      expect_equal(dim(U_hat_norm), c(n_cond, n_voxels))
      for(i in 1:nrow(U_hat_norm)) {
        if (any(!is.na(U_hat_norm[i,])) ) {
          expect_equal(sqrt(sum(U_hat_norm[i,]^2)), 1, tolerance=1e-6)
        }
      }
      
      res_list_avg_folds <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, 
                                                          estimation_method = "average", 
                                                          return_folds = TRUE)
      expect_equal(U_folds_raw, res_list_avg_folds$fold_estimates, tolerance=1e-9)
    }
  )
})

test_that("Original: compute_crossvalidated_means_sl handles single voxel data", {
  n_samples <- 6; n_cond <- 2; n_blocks <- 2; n_voxels <- 1
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  sl_data <- matrix(rnorm(n_samples * n_voxels), nrow = n_samples)
  colnames(sl_data) <- "V1"
  W <- matrix(1,1,1)

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      expect_silent(U_avg <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average"))
      expect_equal(dim(U_avg), c(n_cond, n_voxels))
      expect_silent(U_cn <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "crossnobis", whitening_matrix_W=W))
      expect_equal(dim(U_cn), c(n_cond, n_voxels))
      expect_silent(U_l2 <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "L2_norm"))
      expect_equal(dim(U_l2), c(n_cond, n_voxels))
    }
  )
})

test_that("Original: compute_crossvalidated_means_sl returns NAs if condition missing in training fold", {
  n_samples <- 4; n_cond <- 2; n_voxels <- 3
  # Y_labels: C1, C1, C2, C2
  # custom_block_var: 1, 1, 2, 2 (to make C1 in B1 only, C2 in B2 only)
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks=2, # n_blocks will be updated by custom_block_var
                                  Y_labels = c("C1","C1","C2","C2"),
                                  custom_block_var = c(1,1,2,2))
  cv_spec <- mock_cv_spec_s3(mvpa_des) # .n_folds_val will be 2
  sl_data <- matrix(1:(n_samples*n_voxels), n_samples, n_voxels)

  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = .mock_train_indices,
    .package = "rMVPA",
    {
      res_list <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, "average", return_folds = TRUE)
      
      # Fold 1 trains on B2 data (Y="C2","C2"). C1 is missing.
      expect_true(all(is.na(res_list$fold_estimates["C1", , 1])))
      expect_false(any(is.na(res_list$fold_estimates["C2", , 1])))
      
      # Fold 2 trains on B1 data (Y="C1","C1"). C2 is missing.
      expect_false(any(is.na(res_list$fold_estimates["C1", , 2])))
      expect_true(all(is.na(res_list$fold_estimates["C2", , 2])))
      
      # Final average estimate: C1 (count=1), C2 (count=1). n_folds=2. Both should be NA.
      expect_true(all(is.na(res_list$mean_estimate["C1",])))
      expect_true(all(is.na(res_list$mean_estimate["C2",])))
    }
  )
})

test_that("compute_crossvalidated_means_sl skips empty training folds gracefully", {
  n_samples <- 4; n_cond <- 2; n_blocks <- 2; n_voxels <- 2
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks)
  cv_spec_empty_fold <- mock_cv_spec_s3(mvpa_des)
  
  custom_mock_train_indices <- function(obj, fold_num, ...) {
    if (fold_num == 1) return(integer(0))
    if (inherits(obj, "mock_cv_spec")) return(which(obj$.folds_val != fold_num))
    stop("custom_mock_train_indices error")
  }

  sl_data <- matrix(1:(n_samples*n_voxels), n_samples, n_voxels)
  with_mocked_bindings(
    get_nfolds = .mock_get_nfolds,
    train_indices = custom_mock_train_indices,
    .package = "rMVPA",
    {
      expect_warning(
        res <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec_empty_fold, "average", return_folds = TRUE),
        "Fold 1 has no training samples. Skipping."
      )
      expect_true(all(is.na(res$fold_estimates[,,1])))
      expect_false(any(is.na(res$fold_estimates[,,2])))
    }
  )
}) 