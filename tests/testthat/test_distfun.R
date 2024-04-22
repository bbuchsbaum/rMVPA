library(testthat)

# Test for create_dist Function
test_that("create_dist returns correct structure", {
  dist_obj <- create_dist("euclidean", labels = c("A", "B", "C"))
  expect_true(inherits(dist_obj, "distfun"))
  expect_true(inherits(dist_obj, "euclidean"))
  expect_equal(dist_obj$labels, c("A", "B", "C"))
})

# Tests for pairwise distance functions
test_that("pairwise_dist.correlation calculates correctly", {
  X <- matrix(rnorm(100), 10, 10)
  block <- rep(1:2, each=5)
  dist_obj <- cordist(rep("A", 10), method="pearson")
  result <- pairwise_dist.correlation(dist_obj, X)
  expect_true(is.dist(result))
  expect_equal(length(result), choose(10, 2))
})

test_that("pairwise_dist.euclidean handles empty and single-row matrices", {
  expect_no_error(pairwise_dist.euclidean(list(), matrix(numeric(0), ncol = 10)))
  expect_no_error(pairwise_dist.euclidean(list(), matrix(rnorm(10), nrow = 1)))
})

# Testing Mahalanobis distance with synthetic data
test_that("pairwise_dist.mahalanobis computes correctly", {
  X <- matrix(rnorm(20), 5, 4)
  dist_obj <- list()
  result <- pairwise_dist.mahalanobis(dist_obj, X)
  expect_true(inherits(result, "dist"))
})

test_that("second_order_similarity computes correct scores", {
  # Create a test matrix X with 10 samples and 10 features
  X <- matrix(rnorm(100), 10, 10)
  
  # Create a reference similarity matrix D with the same dimensions as X
  D <- matrix(rnorm(100), 10, 10)
  
  # Define block structure where each block contains exactly two elements
  block <- rep(1:5, each=2)
  
  # Assuming 'create_dist' correctly initializes a Euclidean distance function object
  dist_fun <- create_dist("euclidean", labels = 1:10)
  
  # Compute the similarity scores
  scores <- second_order_similarity(dist_fun, X, D, block, method = "pearson")
  
  # Check the type and length of the output scores
  expect_type(scores, "double")
  expect_length(scores, 10)
  
  # Ensure no scores are NA when blocks are handled correctly
  # This checks that the scores are computed and are not NA due to block handling errors
  expect_false(any(is.na(scores)), "There should be no NA values in the scores if blocks are handled correctly.")
  
  # Test specific behavior for each block
  # Here, we're not expecting any specific outcomes, but you can add these checks based on your function's specifics
  for (i in seq_along(block)) {
    other_block_indices <- which(block != block[i])
    if (length(other_block_indices) > 0) {
      # Test that the score is not NA if there are valid indices to compare
      expect_true(!is.na(scores[i]), sprintf("Score for index %d should not be NA", i))
    }
  }
})



