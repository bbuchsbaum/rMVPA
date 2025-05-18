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

# --- Deterministic dataset: 3 conditions x 2 voxels, two blocks ---
Y_labels <- factor(c("C1","C2","C3","C1","C2","C3"))
block_var <- factor(rep(1:2, each = 3))

mvpa_des <- mock_mvpa_design_cv(n_samples = 6, n_cond = 3, n_blocks = 2, Y_labels = Y_labels)

# Data matrix with condition patterns: C1=c(1,0), C2=c(0,1), C3=c(1,1)
train_mat <- matrix(c(1,0,
                      0,1,
                      1,1,
                      1,0,
                      0,1,
                      1,1), nrow = 6, ncol = 2, byrow = TRUE)

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

# Orthonormal contrast matrix
C_mat <- cbind(
  c(1,0,-1)/sqrt(2),
  c(1,-2,1)/sqrt(6)
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

result_list <- train_model(spec, sl_data = mvpa_dset$train_data, sl_info = sl_info, cv_spec = cv_spec)

# Expected metrics for the first voxel
beta_exp  <- c(C1 = -2.5, C2 = -1.5)
delta_exp <- c(C1 = 0, C2 = sqrt(6)/3)
beta_delta_exp <- beta_exp * delta_exp
beta_delta_norm_exp <- beta_exp * c(0,1)
beta_delta_rel_exp <- beta_delta_exp
composite_exp <- -1.5
recon_exp <- -0.5


test_that("deterministic metrics match manual computation", {
  expect_equal(result_list$beta_only, beta_exp, tolerance = 1e-8)
  expect_equal(result_list$delta_only, delta_exp, tolerance = 1e-8)
  expect_equal(result_list$beta_delta, beta_delta_exp, tolerance = 1e-8)
  expect_equal(result_list$beta_delta_norm, beta_delta_norm_exp, tolerance = 1e-8)
  expect_equal(result_list$beta_delta_reliable, beta_delta_rel_exp, tolerance = 1e-8)
  expect_equal(result_list$composite, setNames(composite_exp, "composite"), tolerance = 1e-8)
  expect_equal(result_list$recon_score, setNames(recon_exp, "recon_score"), tolerance = 1e-8)
})

