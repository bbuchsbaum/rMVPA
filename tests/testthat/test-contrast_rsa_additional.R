context("contrast_rsa additional tests")

library(testthat)
library(rMVPA)

# Helper functions (duplicated from other tests)
mock_mvpa_design_cv <- function(n_samples = 24, n_cond = 4, n_blocks = 3, Y_labels=NULL) {
  if (is.null(Y_labels)) {
    cond_labels <- factor(rep(paste0("Cond", 1:n_cond), each = n_samples / n_cond))
  } else {
    cond_labels <- factor(Y_labels)
    n_cond <- length(levels(cond_labels))
  }
  block_var <- factor(rep(1:n_blocks, length.out = n_samples))
  conditions <- levels(cond_labels)
  design_matrix <- matrix(0, nrow = n_samples, ncol = 1)
  structure(
    list(
      Y = cond_labels,
      block_var = block_var,
      conditions = conditions,
      ncond = n_cond,
      design_matrix = design_matrix,
      samples = 1:n_samples
    ),
    class = c("mvpa_design", "list")
  )
}

mock_cv_spec_s3 <- function(mvpa_design) {
  n_folds_val <- length(unique(mvpa_design$block_var))
  folds_val <- as.integer(mvpa_design$block_var)
  obj <- list(
    .folds_val = folds_val,
    .n_folds_val = n_folds_val
  )
  class(obj) <- c("mock_cv_spec", "cross_validation", "list")
  obj
}

get_nfolds.mock_cv_spec <- function(obj, ...) obj$.n_folds_val
train_indices.mock_cv_spec <- function(obj, fold_num, ...) which(obj$.folds_val != fold_num)

# Define the fixed contrast matrix used in other tests
C_fixed_centered_4x2 <- matrix(c(
  1,  1, -1, -1,  # Contrast 1: (Cond1, Cond2) vs (Cond3, Cond4)
  1, -1,  1, -1   # Contrast 2: (Cond1, Cond3) vs (Cond2, Cond4)
), nrow = 4, ncol = 2, byrow = FALSE)
C_fixed_centered_4x2_scaled <- scale(C_fixed_centered_4x2, center = TRUE, scale = FALSE)
colnames(C_fixed_centered_4x2_scaled) <- c("CFC1", "CFC2")

# -------------------------------------------------------------------------
# Test normalize_delta behaviour
# -------------------------------------------------------------------------

test_that("normalize_delta rescales beta_delta and delta_only", {
  Y_labels <- factor(c("C1","C2","C3","C1","C2","C3"))
  mvpa_des <- mock_mvpa_design_cv(n_samples = 6, n_cond = 3, n_blocks = 2, Y_labels = Y_labels)

  train_mat <- matrix(c(1,0,
                        0,1,
                        1,1,
                        1,0,
                        0,1,
                        1,1), nrow = 6, ncol = 2, byrow = TRUE)

  mvpa_dset <- structure(
    list(
      train_data = train_mat,
      mask = array(1, dim = c(2,2,2)),
      has_test_set = FALSE,
      design = mvpa_des,
      nfeatures = 2
    ),
    class = c("mvpa_dataset", "list")
  )

  C_mat <- cbind(
    c(1,0,-1)/sqrt(2),
    c(1,-2,1)/sqrt(6)
  )
  colnames(C_mat) <- c("C1","C2")
  rownames(C_mat) <- levels(mvpa_des$Y)
  ms_des <- msreve_design(mvpa_des, C_mat)

  spec_raw <- contrast_rsa_model(
    dataset = mvpa_dset,
    design = ms_des,
    output_metric = c("beta_delta","delta_only"),
    normalize_delta = FALSE,
    check_collinearity = FALSE
  )
  cv_spec <- mock_cv_spec_s3(mvpa_des)
  sl_info <- list(center_local_id = 1, center_global_id = 1)
  # Mock function for get_nfolds
  .mock_get_nfolds_contrast <- function(obj, ...) {
    if (inherits(obj, "mock_cv_spec")) {
      return(obj$.n_folds_val)
    }
    stop(".mock_get_nfolds_contrast called with unexpected object type in this test context.")
  }
  
  .mock_train_indices_contrast <- function(obj, fold_num, ...) {
    if (inherits(obj, "mock_cv_spec")) {
      return(which(obj$.folds_val != fold_num))
    }
    stop(".mock_train_indices_contrast called with unexpected object type in this test context.")
  }
  
  res_raw <- with_mocked_bindings(
    get_nfolds = .mock_get_nfolds_contrast,
    train_indices = .mock_train_indices_contrast,
    .package = "rMVPA",
    {
      suppressWarnings(train_model(spec_raw, sl_data = mvpa_dset$train_data, sl_info = sl_info, cv_spec = cv_spec))
    }
  )

  spec_norm <- contrast_rsa_model(
    dataset = mvpa_dset,
    design = ms_des,
    output_metric = c("beta_delta","delta_only"),
    normalize_delta = TRUE,
    check_collinearity = FALSE
  )
  res_norm <- with_mocked_bindings(
    get_nfolds = .mock_get_nfolds_contrast,
    train_indices = .mock_train_indices_contrast,
    .package = "rMVPA",
    {
      suppressWarnings(train_model(spec_norm, sl_data = mvpa_dset$train_data, sl_info = sl_info, cv_spec = cv_spec))
    }
  )

  norm_const <- sqrt(sum(res_raw$delta_only^2))
  expect_equal(res_norm$delta_only, res_raw$delta_only / norm_const, tolerance = 1e-8)
  expect_equal(res_norm$beta_delta, res_raw$beta_delta / norm_const, tolerance = 1e-8)
})

# -------------------------------------------------------------------------
# Test allow_nonorth_composite behaviour
# -------------------------------------------------------------------------

test_that("allow_nonorth_composite returns value with warning", {
  set.seed(987)
  n_samples <- 16; n_cond <- 4; n_blocks <- 2; n_voxels <- 6
  Y_labels <- factor(rep(paste0("Cond", 1:n_cond), each = n_samples / n_cond))
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks, Y_labels)
  base_mat <- matrix(rep(1:n_cond, each = n_samples / n_cond), nrow = n_samples, ncol = n_voxels)
  noise <- matrix(rnorm(n_samples * n_voxels, sd = 0.01), nrow = n_samples)
  data_mat <- base_mat + noise
  dset <- structure(list(train_data = data_mat,
                         mask = array(1, dim = c(2,2,2)),
                         has_test_set = FALSE,
                         design = mvpa_des,
                         nfeatures = n_voxels),
                    class = c("mvpa_dataset", "list"))

  C_nonortho <- C_fixed_centered_4x2_scaled
  rownames(C_nonortho) <- levels(Y_labels)
  ms_des <- msreve_design(mvpa_des, C_nonortho)

  spec <- contrast_rsa_model(dset, ms_des,
                             output_metric = c("composite"),
                             allow_nonorth_composite = TRUE,
                             check_collinearity = FALSE)

  .mock_get_nfolds_contrast2 <- function(obj, ...) {
    if (inherits(obj, "mock_cv_spec")) {
      return(obj$.n_folds_val)
    }
    stop(".mock_get_nfolds_contrast2 called with unexpected object type in this test context.")
  }
  
  .mock_train_indices_contrast2 <- function(obj, fold_num, ...) {
    if (inherits(obj, "mock_cv_spec")) {
      return(which(obj$.folds_val != fold_num))
    }
    stop(".mock_train_indices_contrast2 called with unexpected object type in this test context.")
  }
  
  expect_warning(
    result_list <- with_mocked_bindings(
      get_nfolds = .mock_get_nfolds_contrast2,
      train_indices = .mock_train_indices_contrast2,
      .package = "rMVPA",
      {
        suppressWarnings(train_model(spec,
                                     sl_data = dset$train_data,
                                     sl_info = list(center_local_id = 1, center_global_id = 1, radius = 0, n_voxels = n_voxels),
                                     cv_spec = mock_cv_spec_s3(mvpa_des)))
      }
    ),
    regexp = "not orthonormal"
  )
  expect_false(is.na(result_list$composite))
  expect_null(attr(result_list, "na_reason"))
})

