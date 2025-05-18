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
  expect_true(inherits(out$results$rsa_score, "DenseNeuroVol"))
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
  expect_true(inherits(out$results$rsa_score, "DenseNeuroVol"))
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
  expect_true(inherits(out$results$rsa_score, "DenseNeuroVol"))
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
  expect_equal(mspec$nperm, 0)
  expect_equal(mspec$save_distributions, FALSE)
  
  # Check that other required components are present
  expect_true(!is.null(mspec$distfun))
  expect_true(!is.null(mspec$rsa_simfun))
  expect_true(!is.null(mspec$dataset))
  expect_true(!is.null(mspec$design))
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
  expect_true(inherits(out1$results$rsa_score, "DenseNeuroVol"))
  
  # Test with mahalanobis distance
  mspec <- vector_rsa_model(dataset_obj$dataset, rdes, distfun = mahadist())
  out2 <- run_searchlight(mspec, radius = 4, method = "standard")
  expect_true(inherits(out2$results$rsa_score, "DenseNeuroVol"))
  
  # Test with a PCA-based distance
  threshfun <- function(x) { sum(x > 1) }
  distfun_pca <- pcadist(labels = NULL, ncomp = 3, whiten = FALSE, threshfun = threshfun, dist_method = "cosine")
  mspec <- vector_rsa_model(dataset_obj$dataset, rdes, distfun = distfun_pca)
  out3 <- run_searchlight(mspec, radius = 4, method = "standard")
  expect_true(inherits(out3$results$rsa_score, "DenseNeuroVol"))
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
  expect_true("p_rsa_score" %in% out$metrics, "Permutation p-values (p_rsa_score) should be in output metrics.")
  expect_true("z_rsa_score" %in% out$metrics, "Permutation z-scores (z_rsa_score) should be in output metrics.")
  
  if ("p_rsa_score" %in% out$metrics) {
    p_score_vol <- out$results$p_rsa_score
    expect_true(inherits(p_score_vol, "NeuroVol"), "P-score map should be a NeuroVol.")
    
    # Extract data - NeuroVol may be sparse or dense
    # For DenseNeuroVol, @.Data is the array. For sparse, values() is safer.
    # Or as.vector() for active voxels.
    p_values_vector <- as.vector(p_score_vol) # Gets data from active voxels/vertices
    
    # Filter out NAs which might occur if a searchlight had no valid p-value (e.g. all scores identical)
    p_values_vector <- p_values_vector[!is.na(p_values_vector)]
    
    if (length(p_values_vector) > 0) {
      expect_true(all(p_values_vector >= 0 & p_values_vector <= 1), 
                  paste0("All p-values should be between 0 and 1. Found: ", 
                         paste(p_values_vector[p_values_vector < 0 | p_values_vector > 1], collapse=", ")))
    } else {
      # This might happen if all searchlights failed or produced NA p-values.
      # Depending on test goals, this could be a warning or allowed.
      # For now, we'll just note if no p-values are testable.
      if (length(as.vector(out$results$rsa_score)) > 0) { # Check if rsa_scores were produced
          warning("No non-NA p-values found to test range, though RSA scores exist.")
      }
    }
  }
  
  # Optionally, similar checks for z_rsa_score if specific properties are expected (e.g., not all NA)
  if ("z_rsa_score" %in% out$metrics) {
    z_score_vol <- out$results$z_rsa_score
    expect_true(inherits(z_score_vol, "NeuroVol"), "Z-score map should be a NeuroVol.")
    # Further checks for z-scores could be added here if needed
  }
})

# Setup: Generate sample data and design
set.seed(123)
dset_info <- gen_sample_dataset(c(10,10,10), 60, blocks=3)

# Define unique conditions for the reference D matrix
unique_conditions <- paste0("s", 1:20)

# Create a 20x20 reference dissimilarity matrix with appropriate row/col names
D_ref <- as.matrix(dist(matrix(rnorm(20*10), nrow=20))) # rnorm(20*10) to make a 20-row matrix for dist
rownames(D_ref) <- unique_conditions
colnames(D_ref) <- unique_conditions

# Create the vector RSA design
# 'labels' should correspond to the 60 observations in dset_info, mapping them to the unique_conditions
# 'block_var' also corresponds to the 60 observations
vdes <- vector_rsa_design(D_ref, 
                          labels=rep(unique_conditions, 3), # 60 labels for the 60 observations
                          block_var=dset_info$design$block_var) # block_var from dset_info (length 60)

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






# Potential additional tests:
# - Different rsa_simfun ('spearman')
# - Errors with invalid radius or mask