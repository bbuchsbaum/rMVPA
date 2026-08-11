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

test_that("mvpa_dataset rejects equal-sized images with different spatial geometry", {
  dims <- c(3, 3, 3)
  train_space <- neuroim2::NeuroSpace(
    c(dims, 4), spacing = c(2, 2, 2), origin = c(0, 0, 0)
  )
  shifted_space <- neuroim2::NeuroSpace(
    c(dims, 3), spacing = c(2, 2, 2), origin = c(10, 0, 0)
  )
  mask_space <- neuroim2::NeuroSpace(
    dims, spacing = c(2, 2, 2), origin = c(0, 0, 0)
  )

  train <- neuroim2::NeuroVec(array(rnorm(prod(dims) * 4), c(dims, 4)), train_space)
  test <- neuroim2::NeuroVec(array(rnorm(prod(dims) * 3), c(dims, 3)), shifted_space)
  mask <- neuroim2::NeuroVol(array(1, dims), mask_space)

  expect_error(
    mvpa_dataset(train, test_data = test, mask = mask),
    "Spatial geometry mismatch.*test_data.*origin.*affine"
  )

  aligned_test <- neuroim2::NeuroVec(
    array(rnorm(prod(dims) * 3), c(dims, 3)),
    neuroim2::NeuroSpace(c(dims, 3), spacing = c(2, 2, 2), origin = c(0, 0, 0))
  )
  shifted_mask <- neuroim2::NeuroVol(
    array(1, dims),
    neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2), origin = c(0, 5, 0))
  )
  expect_error(
    mvpa_dataset(train, test_data = aligned_test, mask = shifted_mask),
    "Spatial geometry mismatch.*mask.*origin.*affine"
  )
})

test_that("mvpa_dataset accepts different observation counts in the same spatial geometry", {
  dims <- c(3, 3, 3)
  train_space <- neuroim2::NeuroSpace(
    c(dims, 4), spacing = c(2, 2, 2), origin = c(1, 2, 3)
  )
  test_space <- neuroim2::NeuroSpace(
    c(dims, 2), spacing = c(2, 2, 2), origin = c(1, 2, 3)
  )
  mask_space <- neuroim2::NeuroSpace(
    dims, spacing = c(2, 2, 2), origin = c(1, 2, 3)
  )

  train <- neuroim2::NeuroVec(array(rnorm(prod(dims) * 4), c(dims, 4)), train_space)
  test <- neuroim2::NeuroVec(array(rnorm(prod(dims) * 2), c(dims, 2)), test_space)
  mask <- neuroim2::NeuroVol(array(1, dims), mask_space)

  expect_no_error(mvpa_dataset(train, test_data = test, mask = mask))
})
