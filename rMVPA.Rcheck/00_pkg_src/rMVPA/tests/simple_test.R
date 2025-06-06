# Simple test script
cat("Starting simple test script\n")

# Load required libraries
library(testthat)
library(rMVPA)

# Print a message
cat("Libraries loaded successfully\n")

# Try to run just one test
test_that("Simple test", {
  cat("Running simple test\n")
  expect_true(TRUE)
})

cat("Simple test completed\n")