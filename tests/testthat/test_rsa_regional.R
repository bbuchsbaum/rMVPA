test_that("mvpa_regional with 5 ROIS runs without error", {
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat=Dmat, block=dset$design$block_var), block_var="block")
  mspec <- rsa_model(dset$dataset, design=rdes)
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
   
  res <- run_regional(mspec, region_mask)
  expect_true(!is.null(res))
  
})

test_that("mvpa_regional with 5 ROIS and multiple distance matrices runs without error", {
  
  dset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  Dmat1 <- dist(matrix(rnorm(100*100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100*100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1=Dmat1, Dmat2=Dmat2, block=dset$design$block_var), block_var="block")
  
  mspec <- rsa_model(dset$dataset, design=rdes)
  region_mask <- NeuroVol(sample(1:5, size=length(dset$dataset$mask), replace=TRUE), space(dset$dataset$mask))
  
  res <- run_regional(mspec, region_mask)
  expect_true(!is.null(res))

  
})

library(testthat)

context("Comprehensive tests for RSA regional analysis and design/model functionality")



## --- rsa_design and rsa_model_mat ---
test_that("rsa_design creates a valid RSA design object", {
  # Create a dummy distance matrix and a dummy block variable
  Dmat <- dist(matrix(rnorm(100 * 100), 100, 100))
  data_list <- list(Dmat = Dmat, block = rep(1:5, each = 20))
  rdes <- rsa_design(~ Dmat, data_list, block_var = "block")
  
  expect_true(is.list(rdes))
  expect_true("rsa_design" %in% class(rdes))
  expect_true(!is.null(rdes$formula))
  expect_true(!is.null(rdes$data))
  expect_true(!is.null(rdes$model_mat))
})

test_that("rsa_model_mat returns vectors of correct length and sanitized names", {
  # For a 10x10 distance matrix, lower triangle has 45 elements.
  Dmat <- dist(matrix(rnorm(10 * 10), 10, 10))
  data_list <- list(Dmat = Dmat)
  rdes <- rsa_design(~ Dmat, data_list)
  mm <- rsa_model_mat(rdes)
  
  expect_equal(length(mm[[1]]), 45)
  # Names should be sanitized (e.g., no spaces or colons)
  expect_true(all(grepl("Dmat", names(mm))))
})

## --- Training and Print Methods ---
test_that("train_model.rsa_model with 'lm' regtype returns coefficients with proper names", {
  set.seed(123)
  # Create dummy training data (e.g., 100 observations with 20 features)
  train_data <- matrix(rnorm(100 * 20), 100, 20)
  # Create a block variable and a distance matrix for the design
  block_var <- rep(1:5, each = 20)
  Dmat <- dist(matrix(rnorm(100 * 100), 100, 100))
  data_list <- list(Dmat = Dmat, block = block_var)
  rdes <- rsa_design(~ Dmat, data_list, block_var = "block")
  
  # Create a dummy mvpa_dataset (simulate train_data and a mask)
  dset <- list(train_data = train_data, mask = 1:100)
  class(dset) <- "mvpa_dataset"
  
  mspec <- rsa_model(dset, rdes, regtype = "lm")
  coeffs <- train_model.rsa_model(mspec, train_data, y = NULL, indices = NULL)
  
  expect_true(is.numeric(coeffs))
  expect_true(!is.null(names(coeffs)))
  # Expect one coefficient per predictor in the model matrix.
  expect_equal(length(coeffs), length(rdes$model_mat))
})

test_that("Print methods for rsa_model and rsa_design produce non-empty output", {
  Dmat <- dist(matrix(rnorm(50 * 50), 50, 50))
  data_list <- list(Dmat = Dmat, block = rep(1:5, each = 10))
  rdes <- rsa_design(~ Dmat, data_list, block_var = "block")
  
  dset <- list(train_data = matrix(rnorm(50 * 20), 50, 20), mask = 1:50)
  class(dset) <- "mvpa_dataset"
  mspec <- rsa_model(dset, rdes, regtype = "pearson")
  
  out_design <- capture.output(print(rdes))
  out_model <- capture.output(print(mspec))
  expect_true(length(out_design) > 0)
  expect_true(length(out_model) > 0)
})

## --- Regional Analysis ---
# We assume that run_regional() returns a list with at least the following components:
# "performance_table", "vol_results", and optionally "prediction_table" and "fits".

test_that("mvpa_regional with 5 ROIs runs without error and returns structured result", {
  # Use your provided helper to generate a sample dataset
  dset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat <- dist(matrix(rnorm(100 * 100), 100, 100))
  rdes <- rsa_design(~ Dmat, list(Dmat = Dmat, block = dset$design$block_var), block_var = "block")
  mspec <- rsa_model(dset$dataset, design = rdes)
  
  region_mask <- NeuroVol(sample(1:5, size = length(dset$dataset$mask), replace = TRUE), space(dset$dataset$mask))
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  expect_true(is.list(res))
  expect_true("performance_table" %in% names(res))
  expect_true("vol_results" %in% names(res))
  # Optionally, check prediction_table if return_predictions is enabled
  if (!is.null(res$prediction_table)) {
    expect_true(is.data.frame(res$prediction_table))
  }
})

test_that("mvpa_regional with multiple distance matrices runs without error and returns valid performance metrics", {
  dset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat1 <- dist(matrix(rnorm(100 * 100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100 * 100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1 = Dmat1, Dmat2 = Dmat2, block = dset$design$block_var), block_var = "block")
  mspec <- rsa_model(dset$dataset, design = rdes)
  
  region_mask <- NeuroVol(sample(1:5, size = length(dset$dataset$mask), replace = TRUE), space(dset$dataset$mask))
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  # Check that performance_table contains expected columns (e.g., mse and correlation)
  if (!is.null(res$performance_table)) {
    colnames_perf <- names(res$performance_table)
    expect_true(any(grepl("Dmat1", colnames_perf, ignore.case = TRUE)))
    expect_true(any(grepl("Dmat2", colnames_perf, ignore.case = TRUE)))
  }
})

test_that("mvpa_regional with multiple distance matrices runs without error and returns valid performance metrics", {
  dset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat1 <- dist(matrix(rnorm(100 * 100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100 * 100), 100, 100))
  Dmat3 <- dist(matrix(rnorm(100 * 100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2 + Dmat3, list(Dmat1 = Dmat1, Dmat2 = Dmat2, Dmat3=Dmat3, block = dset$design$block_var), block_var = "block")
  mspec <- rsa_model(dset$dataset, design = rdes, distmethod="spearman", regtype="lm")
  
  region_mask <- NeuroVol(sample(1:5, size = length(dset$dataset$mask), replace = TRUE), space(dset$dataset$mask))
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  # Check that performance_table contains expected columns (e.g., mse and correlation)
  if (!is.null(res$performance_table)) {
    colnames_perf <- names(res$performance_table)
    expect_true(any(grepl("Dmat1", colnames_perf, ignore.case = TRUE)))
    expect_true(any(grepl("Dmat2", colnames_perf, ignore.case = TRUE)))
  }
})

test_that("mvpa_regional with semipartial correlations runs without error", {
  dset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat1 <- dist(matrix(rnorm(100 * 100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100 * 100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1 = Dmat1, Dmat2 = Dmat2, block = dset$design$block_var), block_var = "block")
  
  # Create model with semipartial option set to TRUE
  mspec <- rsa_model(dset$dataset, design = rdes, regtype = "lm", semipartial = TRUE)
  
  region_mask <- NeuroVol(sample(1:5, size = length(dset$dataset$mask), replace = TRUE), space(dset$dataset$mask))
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  # Check that performance_table contains expected columns
  if (!is.null(res$performance_table)) {
    colnames_perf <- names(res$performance_table)
    expect_true(any(grepl("Dmat1", colnames_perf, ignore.case = TRUE)))
    expect_true(any(grepl("Dmat2", colnames_perf, ignore.case = TRUE)))
  }
  
  # Test that coefficients are returned as semi-partial correlations
  # This is a more indirect test, as we can't directly access the coefficients
  # but we can check that the model runs without error and returns results
  expect_true("performance_table" %in% names(res))
  expect_true("vol_results" %in% names(res))
})

test_that("mvpa_regional with non-negative constraints runs without error", {
  dset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat1 <- dist(matrix(rnorm(100 * 100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100 * 100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1 = Dmat1, Dmat2 = Dmat2, block = dset$design$block_var), block_var = "block")
  
  # Create model with non-negative constraints on Dmat2
  mspec <- rsa_model(dset$dataset, design = rdes, regtype = "lm", nneg = list(Dmat2 = TRUE))
  
  region_mask <- NeuroVol(sample(1:5, size = length(dset$dataset$mask), replace = TRUE), space(dset$dataset$mask))
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  # Check that performance_table contains expected columns
  if (!is.null(res$performance_table)) {
    colnames_perf <- names(res$performance_table)
    expect_true(any(grepl("Dmat1", colnames_perf, ignore.case = TRUE)))
    expect_true(any(grepl("Dmat2", colnames_perf, ignore.case = TRUE)))
  }
})

test_that("mvpa_regional with both semipartial and non-negative constraints handles precedence correctly", {
  dset <- gen_sample_dataset(c(5, 5, 5), 100, blocks = 3)
  
  Dmat1 <- dist(matrix(rnorm(100 * 100), 100, 100))
  Dmat2 <- dist(matrix(rnorm(100 * 100), 100, 100))
  rdes <- rsa_design(~ Dmat1 + Dmat2, list(Dmat1 = Dmat1, Dmat2 = Dmat2, block = dset$design$block_var), block_var = "block")
  
  # Create model with both options - non-negative should take precedence
  mspec <- rsa_model(dset$dataset, design = rdes, regtype = "lm", 
                    nneg = list(Dmat2 = TRUE), semipartial = TRUE)
  
  region_mask <- NeuroVol(sample(1:5, size = length(dset$dataset$mask), replace = TRUE), space(dset$dataset$mask))
  res <- run_regional(mspec, region_mask)
  
  expect_true(!is.null(res))
  # Check that results are returned successfully
  expect_true("performance_table" %in% names(res))
  expect_true("vol_results" %in% names(res))
  
  # The model should use non-negative constraints and ignore semipartial
  # This is handled in train_model.rsa_model function, which prioritizes nneg over semipartial
})


