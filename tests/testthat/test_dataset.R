library(testthat)
library(rMVPA)

context("mvpa_dataset validation")

# Test for single voxel dataset

test_that("mvpa_dataset errors for single voxel data", {
  train_array <- array(rnorm(1 * 1 * 1 * 3), dim = c(1, 1, 1, 3))
  class(train_array) <- "NeuroVec"
  mask_array <- array(1, dim = c(1, 1, 1))
  class(mask_array) <- "NeuroVol"
  expect_error(mvpa_dataset(train_array, mask = mask_array),
               "Only 1 voxel detected")
})

# Test for mask with only one active voxel

test_that("mvpa_dataset errors for mask with single active voxel", {
  train_array <- array(rnorm(2 * 1 * 1 * 3), dim = c(2, 1, 1, 3))
  class(train_array) <- "NeuroVec"
  mask_array <- array(c(1, 0), dim = c(2, 1, 1))
  class(mask_array) <- "NeuroVol"
  expect_error(mvpa_dataset(train_array, mask = mask_array),
               "Only 1 active voxel")
})
