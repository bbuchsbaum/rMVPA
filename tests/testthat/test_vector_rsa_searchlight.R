test_that("vector_rsa runs without error and produces valid outputs", {
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
  mspec <- vector_rsa_model(dataset$dataset, rdes, distfun=cordist())
  
  out <- run_searchlight(mspec, radius=4, method="standard")
  expect_true(inherits(out[[1]], "DenseNeuroVol"))
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
  
  out <- run_searchlight.vector_rsa(mspec, radius=4, method="standard")
  expect_true(inherits(out[[1]], "DenseNeuroVol"))
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
  
  out <- run_searchlight.vector_rsa(mspec, radius=4, method="standard")
  expect_true(inherits(out[[1]], "DenseNeuroVol"))
  # Set up parallel processing capabilities
  
})

context("vector_rsa")

## --- Tests for vector_rsa_design ---

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
  # Check that the model spec includes the provided distfun and rsa_simfun
  expect_true(!is.null(mspec$distfun))
  expect_true(mspec$rsa_simfun %in% c("pearson", "spearman"))
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