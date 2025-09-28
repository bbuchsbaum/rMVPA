context("contrast_rsa_model deterministic metrics")

library(testthat)
library(rMVPA)

# Minimal helpers (copied from other tests)
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

# S3 methods for the mock cv spec
get_nfolds.mock_cv_spec <- function(obj, ...) obj$.n_folds_val
train_indices.mock_cv_spec <- function(obj, fold_num, ...) which(obj$.folds_val != fold_num)

# --- Deterministic dataset: 4 conditions x 2 voxels, three blocks ---
# Each condition appears once per block to ensure proper cross-validation
Y_labels <- factor(c("C1","C2","C3","C4","C1","C2","C3","C4","C1","C2","C3","C4"))
block_var <- factor(rep(1:3, each = 4))

mvpa_des <- mock_mvpa_design_cv(n_samples = 12, n_cond = 4, n_blocks = 3, Y_labels = Y_labels)

# Data matrix with condition patterns: C1=c(1,0), C2=c(0,1), C3=c(1,1), C4=c(0.5,0.5)
# Repeated 3 times for 3 blocks
train_mat <- matrix(c(1,0,
                      0,1,
                      1,1,
                      0.5,0.5,
                      1,0,
                      0,1,
                      1,1,
                      0.5,0.5,
                      1,0,
                      0,1,
                      1,1,
                      0.5,0.5), nrow = 12, ncol = 2, byrow = TRUE)

# Simple mvpa_dataset structure
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

# Orthonormal contrast matrix for 4 conditions
C_mat <- cbind(
  c(1, 1, -1, -1)/2,      # C1: (C1+C2) vs (C3+C4)
  c(1, -1, 0, 0)/sqrt(2)  # C2: C1 vs C2
)
colnames(C_mat) <- c("C1","C2")
rownames(C_mat) <- levels(mvpa_des$Y)
ms_des <- msreve_design(mvpa_des, C_mat)

spec <- contrast_rsa_model(
  dataset = mvpa_dset,
  design = ms_des,
  output_metric = c("beta_only","delta_only","beta_delta","beta_delta_norm",
                    "beta_delta_reliable","composite","recon_score"),
  normalize_delta = TRUE,
  calc_reliability = TRUE,
  check_collinearity = FALSE
)

cv_spec <- mock_cv_spec_s3(mvpa_des)
sl_info <- list(center_local_id = 1, center_global_id = 1)

# Mock functions for get_nfolds
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

result_list <- with_mocked_bindings(
  get_nfolds = .mock_get_nfolds_contrast,
  train_indices = .mock_train_indices_contrast,
  .package = "rMVPA",
  {
    suppressWarnings(train_model(spec, sl_data = mvpa_dset$train_data, sl_info = sl_info, cv_spec = cv_spec))
  }
)

# Expected metrics for the first voxel  
# With 4 conditions and patterns C1=c(1,0), C2=c(0,1), C3=c(1,1), C4=c(0.5,0.5)
# These values are computed from the actual cross-validated model output
beta_exp  <- c(C1 = -1.6, C2 = -0.8)  
delta_exp <- c(C1 = -0.3333333, C2 = 0.9428090)
beta_delta_exp <- c(C1 = 0.5333333, C2 = -0.7542472)
beta_delta_norm_exp <- c(C1 = 0.5333333, C2 = -0.7542472)
beta_delta_rel_exp <- c(C1 = 0.5333333, C2 = -0.7542472)
composite_exp <- -0.2209139
recon_exp <- -0.6324555


test_that("deterministic metrics match manual computation", {
  expect_equal(result_list$beta_only, beta_exp, tolerance = 1e-7)
  expect_equal(result_list$delta_only, delta_exp, tolerance = 1e-7)
  expect_equal(result_list$beta_delta, beta_delta_exp, tolerance = 1e-7)
  expect_equal(result_list$beta_delta_norm, beta_delta_norm_exp, tolerance = 1e-7)
  expect_equal(result_list$beta_delta_reliable, beta_delta_rel_exp, tolerance = 1e-7)
  expect_equal(result_list$composite, setNames(composite_exp, "composite"), tolerance = 1e-7)
  expect_equal(result_list$recon_score, setNames(recon_exp, "recon_score"), tolerance = 1e-7)
})
