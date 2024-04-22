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
  
  out <- run_searchlight.vector_rsa(mspec, radius=4, method="standard")
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