library(rMVPA)
library(neuroim2)

test_that("regional feature_rsa_model with direct F matrix runs without error", {
  # Generate a sample dataset: small volume, say 5x5x5, with 100 observations
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a feature matrix F, e.g., 100 observations x 20 features
  Fmat <- matrix(rnorm(100*20), 100, 20)
  labels <- paste0("obs", 1:100)
  
  # Create feature_rsa_design using direct F matrix
  fdes <- feature_rsa_design(F=Fmat, labels=labels, max_comps=3)
  
  # Create a feature_rsa_model, for example using 'pls'
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pls", crossval=blocked_cross_validation(dset$design$block_var))
  
  # Create a region mask with 5 ROIs
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), 
  space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(!is.null(res$performance_table))
})

test_that("regional feature_rsa_model with S-based feature extraction runs without error", {
  # Generate a sample dataset: again 5x5x5 volume, 100 obs
  dset <- gen_sample_dataset(c(6,5,5), 100, blocks=3)
  
  # Create a similarity matrix S: must be symmetric and match number of observations
  obs_features <- matrix(rnorm(100*10), 100, 10)
  S <- tcrossprod(base::scale(obs_features))  # similarity matrix
  labels <- paste0("obs", 1:100)
  
  # Create feature_rsa_design using S
  fdes <- feature_rsa_design(S=S, labels=labels, k=10) # reduce to 5 dims
  
  # Create a feature_rsa_model using pca this time
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pca", crossval=blocked_cross_validation(dset$design$block_var))
  
  # Create a region mask with 5 ROIs
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(!is.null(res$performance_table))
  
  # Check that performance_table has expected columns
  expect_true("mean_correlation" %in% colnames(res$performance_table))
  expect_true("cor_difference" %in% colnames(res$performance_table))
  expect_true("voxel_correlation" %in% colnames(res$performance_table))
  expect_true("mse" %in% colnames(res$performance_table))
  expect_true("r_squared" %in% colnames(res$performance_table))
})

test_that("feature_rsa_model with permutation testing works correctly", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(3,3,3), 50, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(50*10), 50, 10)
  labels <- paste0("obs", 1:50)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model with permutation testing
  mspec <- feature_rsa_model(
    dset$dataset, 
    fdes, 
    method="pca", 
    crossval=blocked_cross_validation(dset$design$block_var),
    nperm=10,  # Small number for testing
    permute_by="observations",
    save_distributions=TRUE
  )
  
  # Create a region mask with just 2 ROIs for faster testing
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  
  # Check that permutation results are included in performance table
  perf_cols <- colnames(res$performance_table)
  expect_true(any(grepl("^p_", perf_cols)))  # p-values
  expect_true(any(grepl("^z_", perf_cols)))  # z-scores
  
  # Check specific permutation columns
  expect_true("p_mean_correlation" %in% perf_cols)
  expect_true("z_mean_correlation" %in% perf_cols)
  expect_true("p_cor_difference" %in% perf_cols)
  expect_true("z_cor_difference" %in% perf_cols)
})

test_that("feature_rsa_model with permute_by='features' works correctly", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(3,3,3), 40, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(40*8), 40, 8)
  labels <- paste0("obs", 1:40)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model with permutation testing by features
  mspec <- feature_rsa_model(
    dset$dataset, 
    fdes, 
    method="pca", 
    crossval=blocked_cross_validation(dset$design$block_var),
    nperm=5,  # Small number for testing
    permute_by="features"  # Permute features instead of observations
  )
  
  # Create a region mask with just 2 ROIs for faster testing
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  
  # Check that permutation results are included in performance table
  perf_cols <- colnames(res$performance_table)
  expect_true(any(grepl("^p_", perf_cols)))
  expect_true(any(grepl("^z_", perf_cols)))
})

test_that("feature_rsa_model with cache_pca=TRUE works correctly", {
  set.seed(123)
  
  # Create a small dataset
  dset <- gen_sample_dataset(c(3,3,3), 30, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(30*6), 30, 6)
  labels <- paste0("obs", 1:30)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model with PCA caching
  mspec <- feature_rsa_model(
    dset$dataset, 
    fdes, 
    method="pca", 
    crossval=blocked_cross_validation(dset$design$block_var),
    cache_pca=TRUE  # Enable PCA caching
  )
  
  # Create a region mask with just 2 ROIs
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(!is.null(res$performance_table))
})

test_that("feature_rsa_model with different max_comps values works correctly", {
  set.seed(123)
  
  # Create a small dataset
  dset <- gen_sample_dataset(c(3,3,3), 40, blocks=2)
  
  # Create a feature matrix with 10 dimensions
  Fmat <- matrix(rnorm(40*10), 40, 10)
  labels <- paste0("obs", 1:40)
  
  # Create feature_rsa_design with different max_comps values
  fdes1 <- feature_rsa_design(F=Fmat, labels=labels, max_comps=3)
  fdes2 <- feature_rsa_design(F=Fmat, labels=labels, max_comps=5)
  
  # Create feature_rsa_models
  mspec1 <- feature_rsa_model(dset$dataset, fdes1, method="pca", 
                             crossval=blocked_cross_validation(dset$design$block_var))
  mspec2 <- feature_rsa_model(dset$dataset, fdes2, method="pca", 
                             crossval=blocked_cross_validation(dset$design$block_var))
  
  # Check that max_comps was properly set
  expect_equal(mspec1$max_comps, 3)
  expect_equal(mspec2$max_comps, 5)
  
  # Create a region mask
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis for both models
  res1 <- run_regional(mspec1, region_mask)
  res2 <- run_regional(mspec2, region_mask)
  
  # Check results
  expect_true(!is.null(res1))
  expect_true(!is.null(res2))
  expect_s3_class(res1, "regional_mvpa_result")
  expect_s3_class(res2, "regional_mvpa_result")
})

test_that("regional feature_rsa_model with pca method runs without error and returns performance", {
  dset <- gen_sample_dataset(c(4,4,4), 80, blocks=4)
  
  # Create an F matrix directly
  Fmat <- matrix(rnorm(80*15), 80, 15)
  labels <- paste0("obs", 1:80)
  
  fdes <- feature_rsa_design(F=Fmat, labels=labels) # no dimension reduction, just use as is
  mspec <- feature_rsa_model(dset$dataset, fdes, method="pca", crossval=blocked_cross_validation(dset$design$block_var))
  
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  # Check if performance metrics are available
  expect_true("performance_table" %in% names(res))
  # Check specific performance metrics
  expect_true("mean_correlation" %in% colnames(res$performance_table))
  expect_true("cor_difference" %in% colnames(res$performance_table))
  expect_true("voxel_correlation" %in% colnames(res$performance_table))
  expect_true("mse" %in% colnames(res$performance_table))
  expect_true("r_squared" %in% colnames(res$performance_table))
})


test_that("can compare feature_rsa with different methods", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(4,4,4), 60, blocks=2)
  
  # Create a feature matrix
  Fmat <- matrix(rnorm(60*10), 60, 10)
  labels <- paste0("obs", 1:60)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels, max_comps=10) # reduce to 5 dims
  
  # Create a region mask with just 2 ROIs for faster testing
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Compare different methods
  methods <- c("pca", "pls")
  results_list <- list()
  
  for (method in methods) {
    print(method)
    # Create model with current method
    mspec <- feature_rsa_model(
      dataset = dset$dataset,
      design = fdes,
      method = method,
      crossval = blocked_cross_validation(dset$design$block_var)
    )
    
    # Run regional analysis
    results_list[[method]] <- run_regional(mspec, region_mask)
    
    # Check results
    expect_true(!is.null(results_list[[method]]))
    expect_s3_class(results_list[[method]], "regional_mvpa_result")
    expect_true(!is.null(results_list[[method]]$performance_table))
  }
  
  # Compare performance metrics across methods
  for (method in methods) {
    # Check that each method has the expected performance metrics
    perf_table <- results_list[[method]]$performance_table
    expect_true("mean_correlation" %in% colnames(perf_table))
    expect_true("cor_difference" %in% colnames(perf_table))
    expect_true("voxel_correlation" %in% colnames(perf_table))
    expect_true("mse" %in% colnames(perf_table))
    expect_true("r_squared" %in% colnames(perf_table))
  }
})

test_that("feature_rsa_model with permutation testing produces valid p-values", {
  set.seed(123)
  
  # Create a small dataset for faster testing
  dset <- gen_sample_dataset(c(3,3,3), 40, blocks=2)
  
  # Create a feature matrix with a strong signal
  X <- matrix(rnorm(40*8), 40, 8)
  B <- matrix(rnorm(8*3), 8, 3)
  Fmat <- X %*% B + matrix(rnorm(40*3, sd=0.1), 40, 3)
  labels <- paste0("obs", 1:40)
  
  # Create feature_rsa_design
  fdes <- feature_rsa_design(F=Fmat, labels=labels)
  
  # Create a feature_rsa_model with permutation testing
  mspec <- feature_rsa_model(
    dset$dataset, 
    fdes, 
    method="pca", 
    crossval=blocked_cross_validation(dset$design$block_var),
    nperm=20,  # Small number for testing
    permute_by="observations"
  )
  
  # Create a region mask with just 1 ROI for faster testing
  region_mask <- NeuroVol(rep(1, length(dset$dataset$mask)), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  expect_s3_class(res, "regional_mvpa_result")
  
  # Check that p-values are between 0 and 1
  p_value_cols <- grep("^p_", colnames(res$performance_table), value=TRUE)
  for (col in p_value_cols) {
    p_values <- res$performance_table[[col]]
    expect_true(all(p_values >= 0 & p_values <= 1))
  }
  
  # Check that z-scores are reasonable
  z_score_cols <- grep("^z_", colnames(res$performance_table), value=TRUE)
  for (col in z_score_cols) {
    z_scores <- res$performance_table[[col]]
    expect_true(all(!is.na(z_scores)))
  }
})