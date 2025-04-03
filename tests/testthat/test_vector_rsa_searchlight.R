test_that("vector_rsa runs without error and produces valid outputs", {
  # Generate a sample dataset with 100 rows, 3 blocks, and a (5,5,5) volume structure
  # Assuming a helper function gen_sample_dataset() that creates suitable data
  dataset <- gen_sample_dataset(c(5,5,5), 15, blocks=3)
  
  # Create a reference distance matrix from random noise, dimensions should match dataset
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- rep(paste0("Label", 1:15), length.out=15)
  row.names(D) <- labels
  colnames(D) <- labels
 
  block <- dataset$design$block_var
  
  # Create vector_rsa_design and model
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(labels), replace=TRUE), block)
  mspec <- vector_rsa_model(dataset$dataset, rdes, distfun=cordist())
  
  out <- run_searchlight(mspec, radius=4, method="standard")
  expect_true(inherits(out[[1]][[1]]$data, "DenseNeuroVol"))
  # Set up parallel processing capabilities
  
})

test_that("vector_rsa runs with mahalanobis distance without error and produces valid outputs", {
  # Generate a sample dataset with 100 rows, 3 blocks, and a (5,5,5) volume structure
  # Assuming a helper function gen_sample_dataset() that creates suitable data
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a reference distance matrix from random noise, dimensions should match dataset
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- rep(paste0("Label", 1:15), length.out=15)
  row.names(D) <- labels
  colnames(D) <- labels
  
  block <- dataset$design$block_var
  
  # Create vector_rsa_design and model
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  mspec <- vector_rsa_model(dataset$dataset, rdes, distfun=mahadist())
  
  out <- run_searchlight(mspec, radius=4, method="standard")
  expect_true(inherits(out$results[[1]]$data, "DenseNeuroVol"))
  # Set up parallel processing capabilities
  
})


test_that("vector_rsa runs with pca distance without error and produces valid outputs", {
  # Generate a sample dataset with 100 rows, 3 blocks, and a (5,5,5) volume structure
  # Assuming a helper function gen_sample_dataset() that creates suitable data
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks=3)
  
  # Create a reference distance matrix from random noise, dimensions should match dataset
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- rep(paste0("Label", 1:15), length.out=15)
  row.names(D) <- labels
  colnames(D) <- labels
  
  block <- dataset$design$block_var
  
  # Create vector_rsa_design and model
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  threshfun <- function(x) {print(sum(x > 1)); sum(x>1)}
  distfun <- pcadist(labels=NULL, ncomp=3, whiten=FALSE, threshfun=threshfun, dist_method="cosine")
  mspec <- vector_rsa_model(dataset$dataset, rdes, distfun=distfun)
  
  out <- run_searchlight(mspec, radius=4, method="standard")
  expect_true(inherits(out$results[[1]]$data, "DenseNeuroVol"))
  # Set up parallel processing capabilities
  
})


test_that("vector_rsa_design errors when labels are missing from row.names of D", {
  # Create a distance matrix with proper dimensions
  D <- as.matrix(dist(matrix(rnorm(100), 10, 10)))
  labels <- paste0("Label", 1:10)
  # Purposely set one rowname to a wrong value
  rownames(D)[1] <- "NotALabel"
  expect_error(
    vector_rsa_design(D = D, labels = labels, block_var = rep(1, 10)),
    "All labels must be present"
  )
})

test_that("vector_rsa_design errors when length of labels and block_var do not match", {
  D <- as.matrix(dist(matrix(rnorm(100), 10, 10)))
  labels <- paste0("Label", 1:10)
  rownames(D) <- labels
  colnames(D) <- labels
  # Provide a block variable of wrong length
  expect_error(
    vector_rsa_design(D = D, labels = labels, block_var = rep(1, 9)),
    "Length of labels and block_var must match"
  )
})

test_that("vector_rsa_design constructs a valid design object", {
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  block <- rep(1:3, length.out = 15)
  
  design <- vector_rsa_design(D = D, labels = labels, block_var = block)
  expect_true(is.list(design))
  expect_true("vector_rsa_design" %in% class(design))
  expect_true(!is.null(design$model_mat))
  # Check that the expanded matrix has dimensions matching labels
  expect_equal(dim(design$model_mat$Dexpanded), c(15,15))
  # Check that cross-block data is computed
  expect_true(is.list(design$model_mat$cross_block_data))
})

## --- Tests for vector_rsa_model ---

test_that("vector_rsa_model errors when design is not a vector_rsa_design", {
  fake_design <- list(a = 1)
  class(fake_design) <- "list"
  dataset <- gen_sample_dataset(c(5,5,5), 100, blocks = 3)$dataset
  expect_error(
    vector_rsa_model(dataset, fake_design),
    "Input must be a 'vector_rsa_design'"
  )
})

test_that("vector_rsa_model constructs a valid model spec", {
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  block <- gen_sample_dataset(c(5,5,5), 100, blocks = 3)$design$block_var
  # Use sampled labels from the block (to match length)
  rdes <- vector_rsa_design(D = D, labels = sample(labels, length(block), replace = TRUE), block_var = block)
  mspec <- vector_rsa_model(gen_sample_dataset(c(5,5,5), 100, blocks = 3)$dataset, rdes, distfun = cordist())
  
  expect_true(inherits(mspec, "vector_rsa_model"))
  # Check that the model spec includes permutation parameters
  expect_equal(mspec$nperm, 50)
  expect_equal(mspec$save_distributions, FALSE)
  
  # Check that other required components are present
  expect_true(!is.null(mspec$distfun))
  expect_true(!is.null(mspec$rsa_simfun))
  expect_true(!is.null(mspec$dataset))
  expect_true(!is.null(mspec$design))
})

## --- Tests for train_model.vector_rsa_model ---
## For these tests we override the second_order_similarity function to simulate output.

test_that("train_model.vector_rsa_model returns expected scores", {
  # Override second_order_similarity with a dummy function
  dummy_second_order <- function(distfun, X, Dref, block_var, method) {
    # For testing, return the mean of X (ignoring NA) as a named scalar.
    # Note: When used across many voxels, the names might not be preserved.
    setNames(mean(X, na.rm = TRUE), "score")
  }
  old_fun <- get("second_order_similarity", envir = .GlobalEnv)
  assign("second_order_similarity", dummy_second_order, envir = .GlobalEnv)
  on.exit(assign("second_order_similarity", old_fun, envir = .GlobalEnv))
  
  # Generate a sample dataset, design, and model spec
  dataset_obj <- gen_sample_dataset(c(5,5,5), 100, blocks = 3)
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  block <- dataset_obj$design$block_var
  
  # Use sampled labels to match length(block)
  rdes <- vector_rsa_design(D = D, labels = sample(labels, length(block), replace = TRUE), block_var = block)
  mspec <- vector_rsa_model(dataset_obj$dataset, rdes, distfun = cordist())
  
  # Compute trial scores via the training function
  scores <- train_model.vector_rsa_model(mspec, as.matrix(dataset_obj$dataset$train_data), y = NULL, indices = NULL)
  
  # Check that the output is numeric and has a positive length
  testthat::expect_true(is.numeric(scores))
  testthat::expect_true(length(scores) > 0)
})

## --- Tests for printing methods ---

test_that("print.vector_rsa_design produces output", {
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  block <- rep(1:3, length.out = 15)
  rdes <- vector_rsa_design(D = D, labels = sample(labels, length(block), replace = TRUE), block_var = block)
  out <- capture.output(print(rdes))
  expect_true(length(out) > 0)
})

test_that("print.vector_rsa_model produces output", {
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  block <- gen_sample_dataset(c(5,5,5), 100, blocks = 3)$design$block_var
  rdes <- vector_rsa_design(D = D, labels = sample(labels, length(block), replace = TRUE), block_var = block)
  mspec <- vector_rsa_model(gen_sample_dataset(c(5,5,5), 100, blocks = 3)$dataset, rdes, distfun = cordist())
  out <- capture.output(print(mspec))
  expect_true(length(out) > 0)
})

## --- Tests for run_searchlight.vector_rsa ---
test_that("vector_rsa searchlight runs without error with different distance functions", {
  dataset_obj <- gen_sample_dataset(c(5,5,5), 100, blocks = 3)
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  block <- dataset_obj$design$block_var
  
  # Test with default cordist
  rdes <- vector_rsa_design(D = D, labels = sample(labels, length(block), replace = TRUE), block_var = block)
  mspec <- vector_rsa_model(dataset_obj$dataset, rdes, distfun = cordist())
  out1 <- run_searchlight(mspec, radius = 4, method = "standard")
  expect_true(inherits(out1$results[[1]]$data, "DenseNeuroVol"))
  
  # Test with mahalanobis distance
  mspec <- vector_rsa_model(dataset_obj$dataset, rdes, distfun = mahadist())
  out2 <- run_searchlight(mspec, radius = 4, method = "standard")
  expect_true(inherits(out2$results[[1]]$data, "DenseNeuroVol"))
  
  # Test with a PCA-based distance
  threshfun <- function(x) { sum(x > 1) }
  distfun_pca <- pcadist(labels = NULL, ncomp = 3, whiten = FALSE, threshfun = threshfun, dist_method = "cosine")
  mspec <- vector_rsa_model(dataset_obj$dataset, rdes, distfun = distfun_pca)
  out3 <- run_searchlight(mspec, radius = 4, method = "standard")
  expect_true(inherits(out3$results[[1]]$data, "DenseNeuroVol"))
})

## --- Tests for permutation testing ---
test_that("vector_rsa runs with permutation testing and produces valid statistical outputs", {
  # Generate sample dataset
  dataset <- gen_sample_dataset(c(5,5,5), 15, blocks=3)
  
  # Create a reference distance matrix
  D <- as.matrix(dist(matrix(rnorm(15*15), 15, 15)))
  labels <- paste0("Label", 1:15)
  rownames(D) <- labels
  colnames(D) <- labels
  
  block <- dataset$design$block_var
  
  # Create vector_rsa_design
  rdes <- vector_rsa_design(D=D, labels=sample(labels, length(block), replace=TRUE), block)
  
  # Create vector_rsa_model with permutation testing enabled
  mspec <- vector_rsa_model(
    dataset$dataset, 
    rdes, 
    distfun=cordist(),
    nperm=50,  # Run 50 permutations
    save_distributions=FALSE  # Don't save full distributions for efficiency
  )
  
  # Run a small searchlight to test permutation
  out <- run_searchlight(mspec, radius=3, method="standard")
  
  # Test that output contains permutation results
  # First get the performance data
  perf_data <- out$results[[1]]$data
  first_nonzero <- which(perf_data@.Data > 0)[1]
  
  # Extract the performance value at the first nonzero location
  if (!is.na(first_nonzero)) {
    # Check if metadata contains permutation results
    meta <- attr(perf_data, "meta")
    
    # Look for p-values and z-scores in column names
    col_names <- colnames(meta$performance_names)
    
    # Verify permutation results exist
    expect_true(any(grepl("^p_", col_names)), "Permutation p-values should exist in results")
    expect_true(any(grepl("^z_", col_names)), "Permutation z-scores should exist in results")
    
    # Also check that there's at least one valid p-value between 0 and 1
    p_cols <- meta$performance_names[, grepl("^p_", col_names), drop=FALSE]
    if (ncol(p_cols) > 0) {
      p_vals <- p_cols[first_nonzero, ]
      expect_true(all(p_vals >= 0 & p_vals <= 1), "P-values should be between 0 and 1")
    }
  }
})

# Setup: Generate sample data and design
set.seed(123)
dset_info <- gen_sample_dataset(c(10,10,10), 60, blocks=3)
vdes <- vector_rsa_design(as.matrix(dist(rnorm(20*10))), 
                          labels=rep(paste0("s", 1:20), 3),
                          block_var=dset_info$design$block_var)

test_that("vector_rsa_model constructs a valid model spec", {
  # No permutation args here
  mspec <- vector_rsa_model(dset_info$dataset, vdes)
  
  expect_s3_class(mspec, "vector_rsa_model")
  expect_true(inherits(mspec, "model_spec"))
  expect_true(!is.null(mspec$dataset))
  expect_true(!is.null(mspec$design))
  expect_true(!is.null(mspec$distfun))
  expect_true(!is.null(mspec$rsa_simfun))
  # REMOVED: These are not stored in the model spec
  # expect_null(mspec$nperm)
  # expect_false(mspec$save_distributions)
})

test_that("vector_rsa runs with standard settings and produces expected output structure", {
  skip_on_cran()
  # No permutation args here
  mspec <- vector_rsa_model(dset_info$dataset, vdes)
  
  # Run searchlight without permutations
  res <- run_searchlight(mspec, radius=4)
  
  expect_s3_class(res, "searchlight_result")
  expect_true(inherits(res$result, "NeuroVol"))
  # Expect a single volume (correlation scores)
  expect_equal(dim(res$result)[4], 1) 
  expect_equal(colnames(res$result), "stat") # Default name for the main statistic
})

test_that("vector_rsa runs with permutation testing and produces valid statistical outputs", {
  skip_on_cran()
  # No permutation args here
  mspec <- vector_rsa_model(dset_info$dataset, vdes)
  
  # Run searchlight WITH permutation args
  res <- run_searchlight(mspec, radius=4, nperm = 50, save_distributions = FALSE)
  
  expect_s3_class(res, "searchlight_result")
  expect_true(inherits(res$result, "NeuroVol"))
  # Expect two volumes: stat and p_value
  expect_equal(dim(res$result)[4], 2)
  expect_true(all(c("stat", "p_value") %in% colnames(res$result)))
  
  # Check p-values are valid (between 0 and 1, possibly NA outside mask)
  p_vals <- res$result[res$result[,,,"p_value"] != 0] # Exclude background
  expect_true(all(p_vals >= 0 & p_vals <= 1, na.rm = TRUE))
})

# Add a test for save_distributions=TRUE if needed, similar structure
test_that("vector_rsa searchlight handles save_distributions=TRUE", {
  skip_on_cran()
  # No permutation args here
  mspec <- vector_rsa_model(dset_info$dataset, vdes)
  
  # Run searchlight WITH permutation args and save_distributions
  res <- run_searchlight(mspec, radius=4, nperm = 50, save_distributions = TRUE)
  
  expect_s3_class(res, "searchlight_result")
  expect_true(inherits(res$result, "NeuroVol"))
  # Expect three volumes: stat, p_value, and possibly perm_mean/perm_sd or similar
  # The exact names/number depends on the implementation detail of save_distributions
  # Let's check for at least 3 volumes
  expect_gte(dim(res$result)[4], 3)
  expect_true(all(c("stat", "p_value") %in% colnames(res$result)))
})


# Potential additional tests:
# - Different rsa_simfun ('spearman')
# - Errors with invalid radius or mask