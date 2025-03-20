library(testthat)

context("Vector RSA regional analysis tests")

test_that("vector_rsa regional analysis with 5 ROIs runs without error", {
  # Generate a sample dataset
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dset$design$block_var
  
  # Create vector_rsa_design and model
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  mspec <- vector_rsa_model(dset$dataset, rdes, distfun=cordist())
  
  # Create a region mask with 5 regions
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Verify that the result is not NULL and has expected components
  expect_true(!is.null(res))
  expect_true(is.list(res))
  expect_true("performance_table" %in% names(res))
  expect_true("vol_results" %in% names(res))
})

test_that("vector_rsa regional analysis works with mahalanobis distance", {
  # Generate a sample dataset
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dset$design$block_var
  
  # Create vector_rsa_design and model with mahalanobis distance
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  mspec <- vector_rsa_model(dset$dataset, rdes, distfun=mahadist())
  
  # Create a region mask with 5 regions
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check that result is not NULL and performance table contains correlation values
  expect_true(!is.null(res))
  if (!is.null(res$performance_table)) {
    # Check if correlation exists as a column in performance_table
    expect_true("correlation" %in% colnames(res$performance_table))
    # Check that correlation values are in a reasonable range (-1 to 1)
    expect_true(all(res$performance_table$correlation >= -1 & res$performance_table$correlation <= 1, na.rm=TRUE))
  }
})

test_that("vector_rsa regional analysis works with PCA-based distance", {
  # Generate a sample dataset
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dset$design$block_var
  
  # Create vector_rsa_design 
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  
  # Create PCA distance function
  threshfun <- function(x) { sum(x > 1) }
  distfun_pca <- pcadist(labels=NULL, ncomp=3, whiten=FALSE, threshfun=threshfun, dist_method="cosine")
  
  # Create model with PCA distance
  mspec <- vector_rsa_model(dset$dataset, rdes, distfun=distfun_pca)
  
  # Create a region mask with 5 regions
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check results
  expect_true(!is.null(res))
  if (!is.null(res$vol_results)) {
    # Check that vol_results contains expected number of volumes
    expect_equal(length(res$vol_results), 1)
  }
})

test_that("vector_rsa regional analysis with permutation testing produces valid statistics", {
  # Generate a sample dataset (smaller for speed)
  dset <- gen_sample_dataset(c(5,5,5), 15, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dset$design$block_var
  
  # Create vector_rsa_design
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  
  # Create model with permutation testing enabled
  mspec <- vector_rsa_model(
    dset$dataset, 
    rdes, 
    distfun=cordist(),
    rsa_simfun="pearson",
    nperm=10,  # Small number for tests to run quickly
    save_distributions=FALSE
  )
  
  # Verify model spec has correct permutation parameters
  expect_equal(mspec$nperm, 10)
  expect_equal(mspec$save_distributions, FALSE)
  
  # Create a region mask with 3 regions (fewer for speed)
  region_mask <- NeuroVol(sample(1:3, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check that result contains permutation statistics
  expect_true(!is.null(res))
  
  if (!is.null(res$performance_table)) {
    # Check for p-value and z-score columns
    col_names <- colnames(res$performance_table)
    expect_true(any(grepl("^p_", col_names)), "Permutation p-values should exist in results")
    expect_true(any(grepl("^z_", col_names)), "Permutation z-scores should exist in results")
    
    # Check that p-values are between 0 and 1
    p_cols <- res$performance_table[, grepl("^p_", col_names), drop=FALSE]
    if (ncol(p_cols) > 0) {
      expect_true(all(p_cols >= 0 & p_cols <= 1, na.rm=TRUE), "P-values should be between 0 and 1")
    }
  }
})

test_that("vector_rsa regional analysis returns correct number of volumes for ROIs", {
  # Generate a sample dataset
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dset$design$block_var
  
  # Create vector_rsa_design and model
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  mspec <- vector_rsa_model(dset$dataset, rdes, distfun=cordist())
  
  # Create a region mask with exactly 4 regions (1,2,3,4)
  mask_data <- sample(c(1,2,3,4), size=length(dset$dataset$mask), replace=TRUE)
  region_mask <- NeuroVol(mask_data, space(dset$dataset$mask))
  
  # Run regional analysis
  res <- run_regional(mspec, region_mask)
  
  # Check that there's one volume per region (4 regions)
  if (!is.null(res$vol_results)) {
    n_regions <- length(unique(mask_data))
    expect_equal(nrow(res$performance_table), n_regions)
  }
})

test_that("vector_rsa regional analysis maintains valid correlation values", {
  # Generate a sample dataset
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dset$design$block_var
  
  # Create vector_rsa_design and model
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  
  # Try both Pearson and Spearman
  mspec_pearson <- vector_rsa_model(dset$dataset, rdes, distfun=cordist(), rsa_simfun="pearson")
  mspec_spearman <- vector_rsa_model(dset$dataset, rdes, distfun=cordist(), rsa_simfun="spearman")
  
  # Create a region mask
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analysis for both models
  res_pearson <- run_regional(mspec_pearson, region_mask)
  res_spearman <- run_regional(mspec_spearman, region_mask)
  
  # Check that correlation values are in valid range (-1 to 1)
  if (!is.null(res_pearson$performance_table)) {
    expect_true(all(res_pearson$performance_table$correlation >= -1 & 
                     res_pearson$performance_table$correlation <= 1, na.rm=TRUE))
  }
  
  if (!is.null(res_spearman$performance_table)) {
    expect_true(all(res_spearman$performance_table$correlation >= -1 & 
                     res_spearman$performance_table$correlation <= 1, na.rm=TRUE))
  }
})

test_that("vector_rsa regional analysis handles save_distributions parameter", {
  # Generate a sample dataset (small for speed)
  dset <- gen_sample_dataset(c(5,5,5), 15, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dset$design$block_var
  
  # Create vector_rsa_design
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  
  # Create two models: one with distributions saved, one without
  mspec_with_dist <- vector_rsa_model(
    dset$dataset, 
    rdes, 
    distfun=cordist(),
    nperm=10,
    save_distributions=TRUE
  )
  
  mspec_without_dist <- vector_rsa_model(
    dset$dataset, 
    rdes, 
    distfun=cordist(),
    nperm=10,
    save_distributions=FALSE
  )
  
  # Create a region mask with 2 regions (for speed)
  region_mask <- NeuroVol(sample(1:2, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  # Run regional analyses
  res_with_dist <- run_regional(mspec_with_dist, region_mask)
  res_without_dist <- run_regional(mspec_without_dist, region_mask)
  
  # Both should run without error and have the same basic structure
  expect_true(!is.null(res_with_dist))
  expect_true(!is.null(res_without_dist))
  
  # Column names should be identical (permutation distributions aren't stored in output table)
  if (!is.null(res_with_dist$performance_table) && !is.null(res_without_dist$performance_table)) {
    expect_equal(colnames(res_with_dist$performance_table), 
                 colnames(res_without_dist$performance_table))
  }
}) 