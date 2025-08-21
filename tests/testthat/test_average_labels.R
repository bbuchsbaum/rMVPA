test_that("average_labels works with basic averaging", {
  library(neuroim2)
  
  # Create test data
  dims <- c(10, 10, 10, 20)  # 20 volumes
  space <- NeuroSpace(dims, c(1,1,1))
  data_array <- array(rnorm(prod(dims)), dims)
  nvec <- NeuroVec(data_array, space)
  
  # Create labels (4 conditions, 5 repetitions each)
  labels <- rep(c("A", "B", "C", "D"), each = 5)
  
  # Test basic averaging
  result <- average_labels(nvec, labels)
  
  # Check dimensions
  expect_equal(dim(result)[1:3], dims[1:3])
  expect_equal(dim(result)[4], 4)  # 4 unique conditions
  
  # Check attributes
  expect_equal(attr(result, "condition_labels"), c("A", "B", "C", "D"))
  expect_equal(as.numeric(attr(result, "n_averaged")), rep(5, 4))
  expect_equal(attr(result, "normalization"), "none")
})

test_that("average_labels works with mask", {
  library(neuroim2)
  
  # Create test data
  dims <- c(10, 10, 10, 12)
  space <- NeuroSpace(dims, c(1,1,1))
  data_array <- array(rnorm(prod(dims)), dims)
  nvec <- NeuroVec(data_array, space)
  
  # Create mask (only center voxels)
  mask_array <- array(FALSE, dims[1:3])
  mask_array[4:7, 4:7, 4:7] <- TRUE
  mask <- NeuroVol(mask_array * 1, NeuroSpace(dims[1:3], c(1,1,1)))
  
  # Create labels
  labels <- rep(c("cond1", "cond2", "cond3"), each = 4)
  
  # Test with mask
  result <- average_labels(nvec, labels, mask = mask)
  
  # Check dimensions
  expect_equal(dim(result)[1:3], dims[1:3])
  expect_equal(dim(result)[4], 3)  # 3 unique conditions
  
  # Check that values outside mask are NA
  # Extract data directly using array indexing
  result_array <- result@.Data
  expect_true(all(is.na(result_array[1,1,1,])))  # Outside mask
  expect_false(all(is.na(result_array[5,5,5,])))  # Inside mask
})

test_that("average_labels normalization methods work", {
  library(neuroim2)
  
  # Create test data with known properties
  dims <- c(5, 5, 5, 10)
  space <- NeuroSpace(dims, c(1,1,1))
  
  # Create data where each volume has different mean and variance
  data_array <- array(0, dims)
  for (i in 1:10) {
    data_array[,,,i] <- rnorm(prod(dims[1:3]), mean = i, sd = i)
  }
  nvec <- NeuroVec(data_array, space)
  
  labels <- rep(c("A", "B"), each = 5)
  
  # Test z-score normalization
  result_z <- average_labels(nvec, labels, normalize = "z", normalize_by = "volume")
  
  # Extract one volume and check it's approximately normalized
  vol1_data <- neuroim2::series(result_z, seq_len(prod(dims[1:3])))[1,]
  vol1_data <- vol1_data[!is.na(vol1_data)]
  expect_true(abs(mean(vol1_data)) < 0.1)  # Mean near 0
  expect_true(abs(sd(vol1_data) - 1) < 0.1)  # SD near 1
  
  # Test unit normalization
  result_unit <- average_labels(nvec, labels, normalize = "unit", normalize_by = "volume")
  
  # Check that each volume has unit norm
  all_data <- neuroim2::series(result_unit, seq_len(prod(dims[1:3])))
  for (i in 1:2) {
    vol_data <- all_data[i,]
    vol_data <- vol_data[!is.na(vol_data)]
    norm_val <- sqrt(sum(vol_data^2))
    expect_true(abs(norm_val - 1) < 0.01)
  }
})

test_that("average_labels handles factor labels correctly", {
  library(neuroim2)
  
  dims <- c(5, 5, 5, 9)
  space <- NeuroSpace(dims, c(1,1,1))
  data_array <- array(rnorm(prod(dims)), dims)
  nvec <- NeuroVec(data_array, space)
  
  # Create factor labels with specific order
  labels <- factor(rep(c("low", "medium", "high"), each = 3),
                   levels = c("low", "medium", "high"))
  
  result <- average_labels(nvec, labels)
  
  # Check that order is preserved
  expect_equal(attr(result, "condition_labels"), c("low", "medium", "high"))
  expect_equal(dim(result)[4], 3)
})

test_that("average_labels return_matrix option works", {
  library(neuroim2)
  
  dims <- c(5, 5, 5, 8)
  space <- NeuroSpace(dims, c(1,1,1))
  data_array <- array(rnorm(prod(dims)), dims)
  nvec <- NeuroVec(data_array, space)
  
  labels <- rep(c("A", "B"), each = 4)
  
  # Get matrix output
  result_mat <- average_labels(nvec, labels, return_matrix = TRUE)
  
  # Check it's a matrix with correct dimensions
  expect_true(is.matrix(result_mat))
  expect_equal(nrow(result_mat), 2)  # 2 conditions
  expect_equal(ncol(result_mat), prod(dims[1:3]))  # Number of voxels
  expect_equal(rownames(result_mat), c("A", "B"))
})

test_that("average_labels handles edge cases", {
  library(neuroim2)
  
  dims <- c(5, 5, 5, 10)
  space <- NeuroSpace(dims, c(1,1,1))
  data_array <- array(rnorm(prod(dims)), dims)
  nvec <- NeuroVec(data_array, space)
  
  # Test with wrong label length
  labels_wrong <- rep("A", 5)  # Only 5 labels for 10 volumes
  expect_error(average_labels(nvec, labels_wrong),
               "Length of labels .* must match number of volumes")
  
  # Test with single condition (all same label)
  labels_single <- rep("A", 10)
  result_single <- average_labels(nvec, labels_single)
  expect_equal(dim(result_single)[4], 1)
  
  # Test with each volume as unique condition (no averaging)
  labels_unique <- as.character(1:10)
  result_unique <- average_labels(nvec, labels_unique)
  expect_equal(dim(result_unique)[4], 10)
})

test_that("normalize_by voxel works correctly", {
  library(neuroim2)
  
  # Create test data
  dims <- c(3, 3, 3, 6)
  space <- NeuroSpace(dims, c(1,1,1))
  
  # Create data where each voxel has different properties
  data_array <- array(rnorm(prod(dims)), dims)
  nvec <- NeuroVec(data_array, space)
  
  labels <- rep(c("A", "B"), each = 3)
  
  # Test voxel-wise z-score normalization
  result <- average_labels(nvec, labels, normalize = "z", normalize_by = "voxel")
  
  # This is harder to test directly, but we can check it runs without error
  expect_equal(dim(result)[4], 2)
  expect_equal(attr(result, "normalization"), "z")
})