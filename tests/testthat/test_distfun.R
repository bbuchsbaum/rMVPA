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
  result <- pairwise_dist(dist_obj, X)
  expect_true(inherits(result, "matrix"))
  expect_equal(length(result), 10*10)
})

test_that("pairwise_dist.euclidean handles empty and single-row matrices", {
  dist_obj <- structure(list(), class="euclidean")
  
  # Empty matrix case
  X <- matrix(numeric(0), ncol = 10)
  result <- pairwise_dist(dist_obj, X)
  expect_equal(dim(result), c(0,0))
  
  # Single row case
  X <- matrix(rnorm(10), nrow = 1)
  result <- pairwise_dist(dist_obj, X)
  expect_equal(dim(result), c(1,1))
  expect_equal(result[1,1], 0)
})

# Testing Mahalanobis distance with synthetic data
test_that("pairwise_dist.mahalanobis computes correctly", {
  dist_obj <- structure(list(), class="mahalanobis")
  X <- matrix(rnorm(20), nrow=4)
  result <- pairwise_dist(dist_obj, X)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(4,4))
  expect_equal(diag(result), rep(0, 4))
})

test_that("pairwise_dist.cordist calculates correctly", {
  dist_obj <- structure(list(method="pearson"), class="cordist")
  
  # Normal case
  X <- matrix(rnorm(20), ncol=4)  # 5x4 matrix
  result <- pairwise_dist(dist_obj, X)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(nrow(X), nrow(X)))  # Should match number of rows
  expect_equal(diag(result), rep(0, nrow(X)))     # Diagonal should match number of rows
  
  # Zero variance case - explicitly create a row with constant values
  X <- rbind(rep(1, 4),                          # Constant row
             matrix(rnorm(12), ncol=4))          # 3 more rows of random data
  expect_warning(
    result <- pairwise_dist(dist_obj, X),
    "the standard deviation is zero"
  )
  expect_true(is.matrix(result))
  expect_true(any(is.na(result)))
})
# Add new tests for pcadist
test_that("pairwise_dist.pcadist works correctly", {
  X <- matrix(rnorm(100), 10, 10)
  
  # Test with default settings
  dist_obj <- pcadist(labels=1:10)
  result1 <- pairwise_dist(dist_obj, X)
  expect_true(inherits(result1, "matrix"))
  expect_equal(dim(result1), c(10,10))
  
  # Test with custom settings
  dist_obj2 <- pcadist(labels=1:10, ncomp=3, whiten=FALSE, dist_method="manhattan")
  result2 <- pairwise_dist(dist_obj2, X)
  expect_true(inherits(result2, "matrix"))
  expect_equal(dim(result2), c(10,10))
})

# Add test for robustmahadist
test_that("pairwise_dist.robustmahadist handles outliers well", {
  # Create data with outliers
  X <- matrix(rnorm(100), 10, 10)
  X[1,] <- X[1,] * 10  # Create an outlier
  
  dist_obj <- robustmahadist()
  result <- pairwise_dist(dist_obj, X)
  
  expect_true(inherits(result, "matrix"))
  expect_equal(dim(result), c(10,10))
  # Distances should be symmetric
  expect_equal(result, t(result))
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

test_that("second_order_similarity handles edge cases", {
  X <- matrix(rnorm(20), 4, 5)
  D <- matrix(c(
    0, 1, 2, 3,
    1, 0, 4, 5,
    2, 4, 0, 6,
    3, 5, 6, 0
  ), nrow=4)
  
  # Test with all samples in same block
  block_same <- rep(1, 4)
  scores_same <- second_order_similarity(eucdist(), X, D, block_same)
  expect_true(all(is.na(scores_same)))
  
  # Test with each sample in different block
  block_diff <- 1:4
  scores_diff <- second_order_similarity(eucdist(), X, D, block_diff)
  expect_false(any(is.na(scores_diff)))
  
  # Test with paired blocks
  block_paired <- rep(1:2, each=2)
  scores_paired <- second_order_similarity(eucdist(), X, D, block_paired)
  expect_equal(sum(is.na(scores_paired)), 0)
})

#' @export
pairwise_dist.default <- function(obj, X, ...) {
  stop(sprintf("pairwise_dist not implemented for objects of class %s", 
               paste(class(obj), collapse=", ")))
}

#' @export
pairwise_dist.cordist <- function(obj, X, ...) {
  if (nrow(X) < 2) {
    return(matrix(0, nrow(X), nrow(X)))
  }
  
  # Check for zero standard deviation
  sds <- apply(X, 2, sd)
  if (any(sds == 0)) {
    warning("the standard deviation is zero")
    # Return matrix with NAs except diagonal
    d <- matrix(NA, nrow(X), nrow(X))
    diag(d) <- 0
    return(d)
  }
  
  # Calculate correlations and convert to distances
  cors <- stats::cor(t(X), method = obj$method)
  1 - cors
}



