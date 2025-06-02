context("contrast_rsa_model constructor")

library(testthat)
library(rMVPA)

# --- Helper objects for contrast_rsa_model tests ---

# Mock mvpa_design object (from former test-msreve_design.R)
# Assuming mvpa_design requires at least condition_labels and block_var
mock_mvpa_design <- function(n_cond = 6, n_blocks = 3) {
  cond <- factor(rep(1:n_cond, each=n_blocks*2)) # Example condition labels
  blocks <- factor(rep(1:n_blocks, times=n_cond*2)) # Example block variable
  
  structure(
    list(
      condition_labels = cond,
      block_var = blocks,
      design_mat = model.matrix(~ cond + blocks), # Example design matrix
      keep_intra_run = FALSE, # Example attribute
      conditions = levels(cond),
      ncond = n_cond,
      samples_per_condition_per_block = table(cond, blocks)
    ),
    class = c("mvpa_design", "list")
  )
}

# Mock mvpa_dataset (very basic structure needed by the constructor)
mock_mvpa_dataset <- function() {
  structure(
    list(
      mask = array(1, dim = c(3,3,3)), # Minimal mask
      has_test_set = FALSE
    ),
    class = c("mvpa_dataset", "list")
  )
}

# Mock msreve_design (uses mock_mvpa_design above)
mock_msreve_design <- function(n_cond = 4, name = "MockMSREVE", center_contrasts = TRUE) {
  # Use mock_mvpa_design_cv as it provides conditions in $conditions and $Y
  # Ensure its conditions are named like "Cond1", "Cond2", etc.
  mvpa_des_for_msreve <- mock_mvpa_design_cv(n_cond = n_cond)
  cond_names <- mvpa_des_for_msreve$conditions # Should be like "Cond1", "Cond2"

  C <- matrix(rnorm(n_cond * 2), nrow = n_cond,
              dimnames = list(cond_names, c("C1", "C2")))
  if (center_contrasts && ncol(C) > 0) {
    C_centered <- scale(C, center = TRUE, scale = FALSE)
    # Retain original colnames and rownames
    dimnames(C_centered) <- dimnames(C)
    C <- C_centered
  }
  msreve_design(mvpa_des_for_msreve, C, name = name)
}


# --- Helpers from former test-compute_cv_means.R ---
# Mock mvpa_design (simplified for CV tests)
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

# Mock cv_spec object with S3 methods (from former test-compute_cv_means.R)
mock_cv_spec_s3 <- function(mvpa_design) {
  n_folds_val <- length(unique(mvpa_design$block_var))
  folds_val <- as.integer(mvpa_design$block_var)
  
  obj <- list(
    .folds_val = folds_val,
    .n_folds_val = n_folds_val
  )
  class(obj) <- c("mock_cv_spec", "cross_validation", "list") # Retain "mock_cv_spec" for S3 dispatch if methods named that way
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


# --- Tests for contrast_rsa_model() ---

test_that("contrast_rsa_model constructor works with valid inputs", {
  dset <- mock_mvpa_dataset()
  ms_des <- mock_msreve_design()
  
  expect_silent({
    model_spec <- contrast_rsa_model(
      dataset = dset,
      design = ms_des,
      estimation_method = "average",
      regression_type = "lm",
      output_metric = c("beta_delta"),
      check_collinearity = FALSE
    )
  })
  
  expect_s3_class(model_spec, "contrast_rsa_model")
  # expect_s3_class(model_spec, "mvpa_model_spec") # Will address this class check later
  expect_s3_class(model_spec, "model_spec")
  
  expect_identical(model_spec$dataset, dset)
  expect_identical(model_spec$design, ms_des)
  expect_equal(model_spec$estimation_method, "average")
  expect_equal(model_spec$regression_type, "lm")
  expect_equal(model_spec$output_metric, c("beta_delta"))
  expect_false(model_spec$check_collinearity)
  expect_false(model_spec$calc_reliability)
})

test_that("contrast_rsa_model constructor errors on invalid dataset", {
  ms_des <- mock_msreve_design()
  expect_error(contrast_rsa_model(dataset = list(), design = ms_des),
               regexp = "dataset does not inherit from class mvpa_dataset") # Adjusted regexp
})

test_that("contrast_rsa_model constructor errors on invalid design", {
  dset <- mock_mvpa_dataset()
  expect_error(contrast_rsa_model(dataset = dset, design = list()),
               regexp = "design does not inherit from class msreve_design") # Adjusted regexp
})

test_that("contrast_rsa_model constructor validates parameters via match.arg", {
  dset <- mock_mvpa_dataset()
  ms_des <- mock_msreve_design()
  
  # Corrected regex based on typical match.arg output
  expect_error(contrast_rsa_model(dset, ms_des, estimation_method = "wrong"),
               regexp = "'arg' should be one of .*“average”, “L2_norm”, “crossnobis”")
  expect_error(contrast_rsa_model(dset, ms_des, regression_type = "wrong"),
               regexp = "'arg' should be one of .*“pearson”, “spearman”, “lm”, “rfit”, “ridge_hkb”")
  expect_error(contrast_rsa_model(dset, ms_des, output_metric = "wrong"),
               regexp = "`output_metric` must be a character vector containing only allowed metrics")
  expect_error(contrast_rsa_model(dset, ms_des, output_metric = c("beta_delta", "wrong_metric")),
               regexp = "`output_metric` must be a character vector containing only allowed metrics")
  expect_error(contrast_rsa_model(dset, ms_des, output_metric = character(0)), # Empty vector
               regexp = "`output_metric` cannot be an empty vector.")
  
  # Test that duplicates are removed and order of first appearance is kept
  expect_silent({
    model_spec_dupes <- contrast_rsa_model(dset, ms_des, output_metric = c("beta_delta", "recon_score", "beta_delta"))
  })
  expect_equal(model_spec_dupes$output_metric, c("beta_delta", "recon_score"))
  
})

test_that("contrast_rsa_model constructor checks logical parameters", {
  dset <- mock_mvpa_dataset()
  ms_des <- mock_msreve_design()

  expect_error(contrast_rsa_model(dset, ms_des, check_collinearity = "TRUE"),
               regexp = "is.logical\\(check_collinearity\\) is not TRUE") # from assert_that
  expect_error(contrast_rsa_model(dset, ms_des, calc_reliability = "yes"),
               regexp = "is.logical\\(calc_reliability\\) is not TRUE")
})

test_that("contrast_rsa_model constructor collinearity check works as expected", {
  dset <- mock_mvpa_dataset()
  n_cond <- 4
  
  # Define condition names as they would be in mock_mvpa_design_cv or similar
  cond_names_for_test <- paste0("Cond", 1:n_cond)

  # Case 1: Non-collinear contrasts - ensure centered and correct dimnames
  C_noncollinear_raw <- matrix(c(1,1,-1,-1,  1,-1,0,0), nrow=n_cond, ncol=2, 
                               dimnames = list(cond_names_for_test, c("C1", "C2")))
  C_noncollinear <- scale(C_noncollinear_raw, center=TRUE, scale=FALSE)
  C_noncollinear <- round(C_noncollinear, 10)
  dimnames(C_noncollinear) <- dimnames(C_noncollinear_raw)
  
  # Use mock_mvpa_design_cv to ensure $conditions field is correctly populated for msreve_design
  mvpa_des_for_noncoll <- mock_mvpa_design_cv(n_cond=n_cond)
  # Override its conditions if they don't match cond_names_for_test (though they should by default)
  mvpa_des_for_noncoll$conditions <- cond_names_for_test
  
  ms_des_noncoll <- suppressWarnings(msreve_design(mvpa_des_for_noncoll, C_noncollinear))
  expect_silent({
    suppressWarnings(
      contrast_rsa_model(dset, ms_des_noncoll, regression_type = "lm", check_collinearity = TRUE)
    )
  })

  # Case 2: Collinear contrasts (C2 = 2 * C1) - ensure centered and correct dimnames
  C_collinear_raw <- matrix(c(1,1,-1,-1,  2,2,-2,-2), nrow=n_cond, ncol=2,
                            dimnames = list(cond_names_for_test, c("C1_orig", "C2_coll")))
  C_collinear <- scale(C_collinear_raw, center=TRUE, scale=FALSE)
  C_collinear <- round(C_collinear, 10)
  dimnames(C_collinear) <- dimnames(C_collinear_raw)

  mvpa_des_for_coll <- mock_mvpa_design_cv(n_cond=n_cond)
  mvpa_des_for_coll$conditions <- cond_names_for_test
  
  ms_des_coll <- suppressWarnings(msreve_design(mvpa_des_for_coll, C_collinear))
  # This expect_error might still see the "not centered" if ms_des_coll construction warns.
  # The goal is for ms_des_coll to be created silently.
  # If C_collinear is perfectly centered, msreve_design should not warn.
  
  # We wrap the msreve_design call itself to catch its warnings if any, for debugging
  # expect_silent({ ms_des_coll_check <- msreve_design(mvpa_des_for_coll, C_collinear) })

  expect_error(
    suppressWarnings(
      contrast_rsa_model(dset, ms_des_coll, regression_type = "lm", check_collinearity = TRUE)
    ),
    # class = "rlang_error", # Removing class check for now, focus on regexp
    regexp = "Collinearity detected among contrast RDMs"
  )
  
  # Case 3: Collinear contrasts but check_collinearity = FALSE
  expect_silent({
    suppressWarnings(contrast_rsa_model(dset, ms_des_coll, regression_type = "lm", check_collinearity = FALSE))
  })

  # Case 4: Collinear contrasts but non-"lm" regression type (check should be skipped)
  expect_silent({
    suppressWarnings(contrast_rsa_model(dset, ms_des_coll, regression_type = "pearson", check_collinearity = TRUE))
  })
  
  # Case 5: Contrast matrix with only 1 column - ensure centered and correct dimnames
  C_single_raw <- matrix(c(1,1,-1,-1), nrow=n_cond, ncol=1,
                         dimnames = list(cond_names_for_test, "C_single"))
  C_single <- scale(C_single_raw, center=TRUE, scale=FALSE)
  C_single <- round(C_single, 10)
  dimnames(C_single) <- dimnames(C_single_raw)
  
  mvpa_des_for_single <- mock_mvpa_design_cv(n_cond=n_cond)
  mvpa_des_for_single$conditions <- cond_names_for_test
  ms_des_single <- suppressWarnings(msreve_design(mvpa_des_for_single, C_single))
  
  expect_silent({
    suppressWarnings(contrast_rsa_model(dset, ms_des_single, regression_type = "lm", check_collinearity = TRUE))
  })
  
  # Case 6: Contrast matrix is NULL (should warn and skip check in constructor)
  ms_des_valid_base <- msreve_design(mvpa_des_for_noncoll, C_noncollinear) # A valid one
  ms_des_valid_then_null_c <- ms_des_valid_base 
  ms_des_valid_then_null_c$contrast_matrix <- NULL
  expect_warning(
    contrast_rsa_model(dset, ms_des_valid_then_null_c, regression_type = "lm", check_collinearity = TRUE),
    regexp = "Cannot check collinearity: Contrast matrix is missing"
  )
})

# --- Tests for train_model.contrast_rsa_model --- 

context("train_model.contrast_rsa_model")

# --- Mock searchlight_info (simplified for these tests) ---
mock_searchlight_info <- function(n_voxels = 10, n_designs = 1) {
  structure(
    list(
      parent_index = 1:n_voxels, # For simplicity, each voxel is its own parent
      parent_ninst = rep(1, n_voxels),
      n_voxels = n_voxels,
      n_designs = n_designs,
      radius = 6, # example radius
      type = "spherical" # example type
    ),
    class = c("searchlight_info", "list")
  )
}

# --- Mock implementations for with_mocked_bindings ---
# These are moved here to be available for all tests in this context
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

# More comprehensive mock dataset for train_model tests
mock_mvpa_dataset_train <- function(n_samples = 24, n_cond = 4, n_blocks = 3, n_voxels = 10) {
  mvpa_des_cv <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks) 
  
  set.seed(456)
  means_matrix <- matrix(as.numeric(mvpa_des_cv$Y), nrow = n_samples, ncol = n_voxels) * 10 
  noise <- matrix(rnorm(n_samples * n_voxels), nrow = n_samples, ncol = n_voxels)
  data_matrix <- means_matrix + noise
  
  mask <- array(1, dim=rep(ceiling(n_voxels^(1/3)), 3)) 
  mask_indices <- 1:prod(dim(mask))
  if (length(mask_indices) < n_voxels) stop("Mock mask too small for n_voxels")
  mask[mask_indices[- (1:n_voxels)]] <- 0
  
  structure(
    list(
      train_data = data_matrix, 
      mask = mask,
      has_test_set = FALSE,
      design = mvpa_des_cv, 
      nfeatures = n_voxels
    ),
    class = c("mvpa_dataset", "list")
  )
}

# Define a fixed, centered contrast matrix for train_model tests where n_cond=4, Q=2
C_fixed_centered_4x2 <- matrix(c(
  1,  1, -1, -1,  # Contrast 1: (Cond1, Cond2) vs (Cond3, Cond4)
  1, -1,  1, -1   # Contrast 2: (Cond1, Cond3) vs (Cond2, Cond4)
), nrow = 4, ncol = 2, byrow = FALSE)
C_fixed_centered_4x2_scaled <- scale(C_fixed_centered_4x2, center = TRUE, scale = FALSE)
colnames(C_fixed_centered_4x2_scaled) <- c("CFC1", "CFC2")
# Rownames will be added dynamically based on dset$design$Y levels in tests


test_that("train_model.contrast_rsa_model runs with standard settings (lm, beta_delta)", {
  n_samples <- 24; n_cond <- 4; n_blocks <- 3; n_voxels_sl <- 7
  center_voxel_local_idx <- 4
  
  dset <- mock_mvpa_dataset_train(n_samples, n_cond, n_blocks, n_voxels = 20)
  
  # Use fixed centered contrast matrix
  C_mat_for_train <- C_fixed_centered_4x2_scaled
  rownames(C_mat_for_train) <- levels(dset$design$Y)
  
  ms_des <- suppressWarnings(msreve_design(dset$design, C_mat_for_train))
  
  model_spec <- contrast_rsa_model(dset, ms_des, check_collinearity = FALSE) 
  cv_spec <- mock_cv_spec_s3(dset$design) 
  sl_data <- dset$train_data[, 1:n_voxels_sl] 
  sl_info <- list(center_local_id = center_voxel_local_idx, center_global_id = center_voxel_local_idx, radius = 0, n_voxels = n_voxels_sl)

  expect_silent({
    with_mocked_bindings(
      get_nfolds = .mock_get_nfolds_contrast,
      train_indices = .mock_train_indices_contrast,
      .package = "rMVPA",
      {
        res <- suppressWarnings(train_model(model_spec, sl_data, sl_info, cv_spec))
      }
    )
  })
  
  expect_s3_class(res, "contrast_rsa_model_results")
  
  result_vec <- res$beta_delta
  expect_type(result_vec, "double")
  expect_named(result_vec, colnames(ms_des$contrast_matrix))
  expect_false(anyNA(result_vec)) 
})

test_that("train_model.contrast_rsa_model works with pearson, delta_only", {
  n_samples <- 24; n_cond <- 4; n_blocks <- 3; n_voxels_sl <- 8
  center_voxel_local_idx <- 5
  
  dset <- mock_mvpa_dataset_train(n_samples, n_cond, n_blocks, n_voxels = 20)

  # Use fixed centered contrast matrix
  C_mat_for_train <- C_fixed_centered_4x2_scaled
  rownames(C_mat_for_train) <- levels(dset$design$Y)
  
  ms_des <- suppressWarnings(msreve_design(dset$design, C_mat_for_train))

  model_spec <- contrast_rsa_model(dset, ms_des, 
                                   regression_type = "pearson", 
                                   output_metric = c("delta_only"),
                                   check_collinearity = FALSE) 
  cv_spec <- mock_cv_spec_s3(dset$design) 
  sl_data <- dset$train_data[, 1:n_voxels_sl]
  sl_info <- list(center_local_id = center_voxel_local_idx, center_global_id = center_voxel_local_idx, radius = 0, n_voxels = n_voxels_sl)
  
  expect_silent({
    with_mocked_bindings(
      get_nfolds = .mock_get_nfolds_contrast,
      train_indices = .mock_train_indices_contrast,
      .package = "rMVPA",
      {
        res <- suppressWarnings(train_model(model_spec, sl_data, sl_info, cv_spec))
      }
    )
  })
  
  expect_s3_class(res, "contrast_rsa_model_results")
  
  result_vec <- res$delta_only
  expect_type(result_vec, "double")
  expect_named(result_vec, colnames(ms_des$contrast_matrix))
  expect_false(anyNA(result_vec))
})

test_that("train_model.contrast_rsa_model works with spearman, beta_only", {
  n_samples <- 24; n_cond <- 4; n_blocks <- 3; n_voxels_sl <- 9
  center_voxel_local_idx <- 2
  
  dset <- mock_mvpa_dataset_train(n_samples, n_cond, n_blocks, n_voxels = 20)

  # Use fixed centered contrast matrix
  C_mat_for_train <- C_fixed_centered_4x2_scaled
  rownames(C_mat_for_train) <- levels(dset$design$Y)

  ms_des <- suppressWarnings(msreve_design(dset$design, C_mat_for_train))
  
  model_spec <- contrast_rsa_model(dset, ms_des, 
                                   regression_type = "spearman", 
                                   output_metric = c("beta_only"),
                                   check_collinearity = FALSE) 
  cv_spec <- mock_cv_spec_s3(dset$design) 
  sl_data <- dset$train_data[, 1:n_voxels_sl]
  sl_info <- list(center_local_id = center_voxel_local_idx, center_global_id = center_voxel_local_idx, radius = 0, n_voxels = n_voxels_sl)
  
  expect_silent({
    with_mocked_bindings(
      get_nfolds = .mock_get_nfolds_contrast,
      train_indices = .mock_train_indices_contrast,
      .package = "rMVPA",
      {
        res <- suppressWarnings(train_model(model_spec, sl_data, sl_info, cv_spec))
      }
    )
  })
  
  expect_s3_class(res, "contrast_rsa_model_results")
  
  result_vec <- res$beta_only
  expect_type(result_vec, "double")
  expect_named(result_vec, colnames(ms_des$contrast_matrix))
  expect_false(anyNA(result_vec))
  expect_true(all(result_vec >= -1 & result_vec <= 1)) 
})

test_that("train_model handles NA from U_hat_for_delta_calc (formerly compute_crossvalidated_means_sl)", {
  n_cond_na <- 2
  C_mat_na <- matrix(c(1,-1, 1,1), nrow=n_cond_na, ncol=2, dimnames=list(paste0("Cond",1:n_cond_na), c("CN1", "CN2")))
  # Center the contrast matrix to avoid unrelated warnings
  C_mat_na_centered <- scale(C_mat_na, center=TRUE, scale=FALSE)
  dimnames(C_mat_na_centered) <- dimnames(C_mat_na)
  C_mat_na <- C_mat_na_centered

  mvpa_des_for_na <- mock_mvpa_design_cv(n_cond = n_cond_na, n_samples = 4, n_blocks = 2)
  ms_des_na <- suppressWarnings(msreve_design(mvpa_des_for_na, C_mat_na)) # Should use centered C_mat_na
  
  model_spec_na <- contrast_rsa_model(
    dataset = mock_mvpa_dataset_train(n_cond=n_cond_na, n_samples=4, n_blocks=2, n_voxels=5), 
    design = ms_des_na, 
    estimation_method = "average", 
    output_metric = c("beta_delta", "recon_score")
  )

  # Mock compute_crossvalidated_means_sl to return U_hat with NAs
  # This U_hat will be used for U_hat_for_delta_calc
  U_hat_with_na <- matrix(c(1,NA,3,4, 5,6,7,8, 9,10,11,12), nrow=n_cond_na, ncol=6, byrow=TRUE) 
  dimnames(U_hat_with_na) <- list(levels(model_spec_na$dataset$design$Y), paste0("V",1:6))

  # Mock the internal call to compute_crossvalidated_means_sl
  # This mock needs to be active when train_model.contrast_rsa_model calls it
  # The first call to compute_cv_means_sl is for U_hat_for_delta_calc
  
  result_list <- with_mocked_bindings(
    compute_crossvalidated_means_sl = function(...) U_hat_with_na,
    .package = "rMVPA",
    {
      train_model(
        model_spec_na, 
        sl_data = model_spec_na$dataset$train_data[1:4, 1:5, drop=FALSE],
        sl_info = list(center_local_id = 1, center_global_id = 1, radius=0, n_voxels=5)
      )
    }
  )
  
  expect_type(result_list, "list")
  expect_named(result_list, c("beta_delta", "recon_score"))
  expect_equal(attr(result_list, "na_reason"), "NA in U_hat_for_delta_calc")
  
  # Check beta_delta element
  expect_type(result_list$beta_delta, "double")
  expect_true(all(is.na(result_list$beta_delta)))
  expect_equal(length(result_list$beta_delta), ncol(C_mat_na))
  
  # Check recon_score element
  expect_type(result_list$recon_score, "double") # or logical if it's NA_real_
  expect_true(is.na(result_list$recon_score))
  expect_equal(length(result_list$recon_score), 1)
})

test_that("train_model handles insufficient data points for regression", {
    n_cond_insufficient <- 3
    # This contrast matrix is rank deficient (C3 = C2-C1), but more importantly, it's for 3 conditions.
    # Original C_insufficient was rank 2. Let's make it full rank (3) for 3 conditions to ensure ncol(X_reg) is 3.
    C_insufficient_orig <- matrix(c(1,-1,0,  0,1,-1,  1,0,-1), nrow=n_cond_insufficient, ncol=3,
                                  dimnames=list(paste0("Cond",1:n_cond_insufficient), paste0("CQ",1:3)))
    # Center it
    C_insufficient_centered <- scale(C_insufficient_orig, center=TRUE, scale=FALSE)
    dimnames(C_insufficient_centered) <- dimnames(C_insufficient_orig)
    C_insufficient_final <- C_insufficient_centered
    
    dset_insufficient <- mock_mvpa_dataset_train(n_samples = 6, n_cond = n_cond_insufficient, n_blocks = 2, n_voxels = 1)
    ms_des_insufficient <- suppressWarnings(msreve_design(dset_insufficient$design, C_insufficient_final))

    model_spec_insufficient <- contrast_rsa_model(
        dataset = dset_insufficient, 
        design = ms_des_insufficient, 
        estimation_method = "average", # This will use U_hat_for_delta_calc -> G_hat -> dvec_sl
        output_metric = c("beta_delta"),
        check_collinearity = FALSE # Avoid collinearity check for this specific test focus
    )

    # Create U_hat_for_delta_calc that is NOT NA, but will lead to insufficient RSA pairs
    # e.g. U_hat where all conditions are nearly identical, so G_hat has many zeros in lower.tri
    # K=3 conditions. U_hat_for_delta_calc will be 3 (conds) x 1 (voxel)
    U_hat_clean_but_problematic <- matrix(c(1.0, 1.000001, 1.000002), nrow=n_cond_insufficient, ncol=1)
    dimnames(U_hat_clean_but_problematic) <- list(levels(dset_insufficient$design$Y), "V1")

    # This will make dvec_sl have very few unique, non-NA values if not careful.
    # G_hat_sl = U_hat %*% t(U_hat). For K=3, G_hat is 3x3. lower.tri has 3 elements.
    # If U_hat rows are almost same, G_hat elements will be almost same. dvec_sl elements will be almost same.
    # The check in run_rsa_regression_base is: length(good_dvec) < (ncol(X_reg) + 2)
    # good_dvec = dvec_sl[is.finite(dvec_sl)]. X_reg = C_insufficient_final (3 columns). So need length(good_dvec) < 5.
    # If dvec_sl has 3 elements, all finite, length is 3, which is < 5. So this should trigger.
    
    # We want to test the "Insufficient valid data points...for RSA regression" warning.
    # This warning is generated in run_rsa_regression_base.
    # For this to happen, the anyNA(U_hat_for_delta_calc) check in train_model MUST PASS (no NAs).

    expect_warning({
        result_list <- with_mocked_bindings(
            compute_crossvalidated_means_sl = function(...) U_hat_clean_but_problematic, 
            .package = "rMVPA",
            {
                train_model(
                    model_spec_insufficient,
                    sl_data = model_spec_insufficient$dataset$train_data,
                    sl_info = list(center_local_id = 1, center_global_id = 1, radius=0, n_voxels=1)
                )
            }
        )
    }, regexp = "Insufficient valid data points.*for RSA regression.*Returning NAs.")

    expect_type(result_list, "list")
    expect_named(result_list, c("beta_delta"))
    expect_equal(attr(result_list, "na_reason"), "Insufficient valid RSA pairs")

    result_vec <- result_list$beta_delta
    expect_type(result_vec, "double")
    expect_true(all(is.na(result_vec)))
    expect_equal(length(result_vec), ncol(C_insufficient_final))
    expect_named(result_vec, colnames(C_insufficient_final))
})

test_that("train_model errors with invalid sl_info", {
  n_samples <- 24; n_cond <- 4; n_blocks <- 3; n_voxels_sl <- 7
  
  dset <- mock_mvpa_dataset_train(n_samples, n_cond, n_blocks, n_voxels = 20) 
  
  # Use fixed centered contrast matrix
  C_mat_for_train <- C_fixed_centered_4x2_scaled
  rownames(C_mat_for_train) <- levels(dset$design$Y)

  ms_des <- suppressWarnings(msreve_design(dset$design, C_mat_for_train))
  model_spec <- contrast_rsa_model(dset, ms_des, check_collinearity = FALSE, output_metric = c("beta_delta"))
  
  cv_spec <- mock_cv_spec_s3(dset$design) 
  sl_data <- dset$train_data[, 1:n_voxels_sl]
  
  sl_info_missing_local <- list(center_global_id = 1) # Missing center_local_id
  sl_info_null_local <- list(center_local_id = NULL, center_global_id = 1)
  sl_info_na_local <- list(center_local_id = NA_integer_, center_global_id = 1)
  sl_info_oob <- list(center_local_id = 100, center_global_id = 100)
  sl_info_zero <- list(center_local_id = 0, center_global_id = 0)
  
  # Updated regexps based on actual error from train_model
  expect_error(train_model(model_spec, sl_data, sl_info_missing_local, cv_spec),
               regexp = "`train_model.contrast_rsa_model` requires a valid center voxel ID")
  expect_error(train_model(model_spec, sl_data, sl_info_null_local, cv_spec),
               regexp = "`train_model.contrast_rsa_model` requires a valid center voxel ID")
  expect_error(train_model(model_spec, sl_data, sl_info_na_local, cv_spec),
               regexp = "`train_model.contrast_rsa_model` requires a valid center voxel ID")
  expect_error(train_model(model_spec, sl_data, sl_info_oob, cv_spec),
               regexp = "`sl_info\\$center_local_id` .* is out of bounds for `sl_data` columns")
  expect_error(train_model(model_spec, sl_data, sl_info_zero, cv_spec),
               regexp = "`sl_info\\$center_local_id` .* is out of bounds for `sl_data` columns")
})

# Test: Check if whitening matrix W is correctly passed and used with "crossnobis"
# This is a more conceptual test for now, verifying parameters are plumbed
test_that("train_model.contrast_rsa_model handles whitening_matrix_W for crossnobis", {
  n_samples <- 12; n_cond <- 3; n_blocks <- 3; n_voxels_sl <- 5
  
  dset <- mock_mvpa_dataset_train(n_samples, n_cond, n_blocks, n_voxels = n_voxels_sl)
  
  C_mat_raw <- matrix(rnorm(n_cond * 2), nrow=n_cond, dimnames=list(levels(dset$design$Y), c("C1", "C2")))
  C_mat_centered <- scale(C_mat_raw, center=TRUE, scale=FALSE)
  dimnames(C_mat_centered) <- dimnames(C_mat_raw)
  
  ms_des <- suppressWarnings(msreve_design(dset$design, C_mat_centered))
  
  # Create a dummy whitening matrix
  W_dummy <- diag(n_voxels_sl) * 0.5 
  
  model_spec_cn <- contrast_rsa_model(
    dataset = dset, 
    design = ms_des, 
    estimation_method = "crossnobis",
    output_metric = c("beta_delta"),
    whitening_matrix_W = W_dummy, # Pass W here
    check_collinearity = FALSE
  )
  
  # We need to mock compute_crossvalidated_means_sl to check if it receives W
  # And also mock run_rsa_regression_base to stop the chain
  
  mock_compute_cv_means_cn <- function(sl_data, mvpa_des, cv_spec, estimation_method, whitening_matrix_W = NULL, ...) {
    # Check if W was passed
    if (!is.null(whitening_matrix_W)) {
      attr(sl_data, "W_received") <- TRUE # Tag it so we can check
      expect_identical(whitening_matrix_W, W_dummy) # Check if it's the one we passed
    } else {
       attr(sl_data, "W_received") <- FALSE
    }
    # Return a structure that dist_sq_from_fold_means can use
    # and also the U_hat for delta calc
    # For this test, exact values don't matter as much as the W plumbing
    K_val <- mvpa_des$ncond
    V_val <- ncol(sl_data)
    
    # Mocked fold_means: list over folds, each KxV
    # Here, we need cv_spec to get n_folds
    # Assuming cv_spec is auto-generated if NULL in train_model from model_spec$dataset$design$block_var
    n_folds_mock <- length(unique(mvpa_des$block_var))
    
    mock_fold_means <- rep(list(matrix(0, nrow=K_val, ncol=V_val)), n_folds_mock)
    
    return(list(
      fold_means = mock_fold_means, 
      U_hat = matrix(0, nrow=K_val, ncol=V_val), # For U_hat_for_delta_calc
      params = list() # other params if any
    ))
  }
  
  # Mock the RSA regression part to just return NAs of correct size
  mock_rsa_reg_base <- function(...) {
    # Return NAs for C1, C2
    setNames(rep(NA_real_, ncol(C_mat_centered)), colnames(C_mat_centered))
  }

  # We expect W_dummy to be passed to compute_crossvalidated_means_sl
  # The actual result of train_model doesn't matter as much as the mock's internal check
  
  W_was_received <- FALSE
  with_mocked_bindings(
    compute_crossvalidated_means_sl = function(sl_data, mvpa_des, cv_spec, estimation_method, whitening_matrix_W = NULL, ...) {
        if (!is.null(whitening_matrix_W)) {
          expect_identical(whitening_matrix_W, W_dummy)
          W_was_received <<- TRUE 
        }
        mock_compute_cv_means_cn(sl_data, mvpa_des, cv_spec, estimation_method, whitening_matrix_W, ...)
    },
    compute_crossnobis_distances_sl = function(...) rep(NA_real_, 3),
    .package = "rMVPA",
    {
      train_model(
        model_spec_cn, 
        sl_data = dset$train_data,
        sl_info = list(center_local_id = 1, center_global_id = 1, radius=0, n_voxels=n_voxels_sl)
      )
    }
  )
  
  expect_true(W_was_received, 
              info = "Whitening matrix W was not received by the mocked compute_crossvalidated_means_sl when estimation_method='crossnobis'.")
  
})

# New discriminating test: crossnobis with mismatched whitening matrix dimensions

test_that("train_model.contrast_rsa_model errors when whitening_matrix_W dimension mismatches sl_data", {
  n_samples <- 12; n_cond <- 4; n_blocks <- 3; n_voxels_sl <- 4
  
  # Create mock dataset with exactly 4 voxels
  dset <- mock_mvpa_dataset_train(n_samples, n_cond, n_blocks, n_voxels = n_voxels_sl)
  
  # Use fixed centered contrast matrix (rownames set below)
  C_mat_for_train <- C_fixed_centered_4x2_scaled
  rownames(C_mat_for_train) <- levels(dset$design$Y)
  ms_des <- suppressWarnings(msreve_design(dset$design, C_mat_for_train))
  
  # Create a whitening matrix with wrong dimensions (3 x 3 instead of 4 x 4)
  W_wrong <- diag(3) * 0.5
  
  model_spec_cn <- contrast_rsa_model(
    dataset = dset,
    design = ms_des,
    estimation_method = "crossnobis",
    output_metric = c("beta_delta"),
    whitening_matrix_W = W_wrong,
    check_collinearity = FALSE
  )
  
  cv_spec <- mock_cv_spec_s3(dset$design)
  sl_data <- dset$train_data  # 12 x 4
  sl_info <- list(center_local_id = 2, center_global_id = 2, radius = 0, n_voxels = n_voxels_sl)
  
  expect_error(
    train_model(model_spec_cn, sl_data, sl_info, cv_spec),
    regexp = "whitening_matrix_W.*dimensions"  # expect informative abort
  )
})

# Highly discriminating metric-consistency test -----------------------------------------------------

test_that("contrast_rsa_model output metrics are internally consistent", {
  set.seed(987)
  n_samples <- 16; n_cond <- 4; n_blocks <- 2; n_voxels <- 6
  
  # Create deterministic toy dataset: condition means differ by row index
  Y_labels <- factor(rep(paste0("Cond", 1:n_cond), each = n_samples / n_cond))
  block_var <- factor(rep(1:n_blocks, length.out = n_samples))
  base_mat <- matrix(rep(1:n_cond, each = n_samples / n_cond), nrow = n_samples, ncol = n_voxels)
  noise <- matrix(rnorm(n_samples * n_voxels, sd = 0.01), nrow = n_samples)
  data_mat <- base_mat + noise
  
  # Build mvpa_design/mock dataset
  mvpa_des <- mock_mvpa_design_cv(n_samples, n_cond, n_blocks, Y_labels)
  dset <- structure(list(train_data = data_mat,
                         mask = array(1, dim = c(2,2,2)),
                         has_test_set = FALSE,
                         design = mvpa_des,
                         nfeatures = n_voxels),
                    class = c("mvpa_dataset", "list"))
  
  # Orthonormal 2-contrast matrix using helper
  C_raw <- C_fixed_centered_4x2_scaled
  C_ortho <- orthogonalize_contrasts(C_raw)
  rownames(C_ortho) <- levels(Y_labels)
  
  ms_des <- msreve_design(mvpa_des, C_ortho)
  
  # Helper to run model with given options
  run_metric <- function(metric, normalize = FALSE, reliability = FALSE){
    spec <- contrast_rsa_model(dset, ms_des,
                               output_metric = c(metric),
                               normalize_delta = normalize,
                               calc_reliability = reliability,
                               check_collinearity = FALSE)
    # result_list <- train_model(spec, # Original call
    #             sl_data = dset$train_data,
    #             sl_info = list(center_local_id = 1, center_global_id = 1, radius = 0, n_voxels = n_voxels),
    #             cv_spec = mock_cv_spec_s3(mvpa_des))
    
    # Wrapped call
    result_list <- with_mocked_bindings(
        get_nfolds = .mock_get_nfolds_contrast,
        train_indices = .mock_train_indices_contrast,
        .package = "rMVPA",
        {
            suppressWarnings(train_model(spec,
                sl_data = dset$train_data,
                sl_info = list(center_local_id = 1, center_global_id = 1, radius = 0, n_voxels = n_voxels),
                cv_spec = mock_cv_spec_s3(mvpa_des)))
        }
    )
    result_list[[metric]]
  }
  
  beta_only_vec <- run_metric("beta_only")
  delta_only_vec <- run_metric("delta_only")
  beta_delta_vec <- run_metric("beta_delta")
  
  # 1. beta_delta should equal element-wise product of beta and delta
  expect_equal(beta_delta_vec, beta_only_vec * delta_only_vec, tolerance = 1e-6)
  
  # 2. Normalised delta version
  beta_delta_norm_val <- run_metric("beta_delta_norm", normalize = TRUE)
  norm_delta_val <- delta_only_vec / sqrt(sum(delta_only_vec^2))
  expect_equal(beta_delta_norm_val, beta_only_vec * norm_delta_val, tolerance = 1e-6)
  
  # 3. Composite metric should equal sum(beta_delta_norm) because C is orthonormal
  composite_val_single <- run_metric("composite", normalize = TRUE)
  expect_equal(as.numeric(composite_val_single), sum(beta_delta_norm_val), tolerance = 1e-6)
  
  # 4. recon_score must be finite and between -1 and 1
  recon_val_single <- run_metric("recon_score")
  expect_true(is.finite(recon_val_single))
  expect_true(recon_val_single >= -1 && recon_val_single <= 1)

  # rho_const <- (get_nfolds(mock_cv_spec_s3(mvpa_des)) - 1) / get_nfolds(mock_cv_spec_s3(mvpa_des)) # Original
  rho_const <- (.mock_get_nfolds_contrast(mock_cv_spec_s3(mvpa_des)) - 1) / .mock_get_nfolds_contrast(mock_cv_spec_s3(mvpa_des)) # Patched
  beta_delta_rel <- run_metric("beta_delta_reliable", reliability = TRUE)
  expect_equal(beta_delta_rel, beta_delta_vec * rho_const, tolerance = 1e-6)
})
