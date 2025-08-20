test_that("searchlight handles filtered center voxel gracefully", {
  library(neuroim2)
  
  # Create test data using pattern from test_mvpa_searchlight.R
  D <- c(5, 5, 5)  # Small 3D volume
  nobs <- 10  # Number of observations
  nlevels <- 2  # Binary classification
  
  # Create 4D data array with a center voxel that has zero variance
  mat <- array(rnorm(prod(D) * nobs), c(D, nobs))
  
  # Set center voxel (position 3,3,3 in 5x5x5) to constant value
  center_pos <- c(3, 3, 3)
  mat[center_pos[1], center_pos[2], center_pos[3], ] <- 5  # All same = zero variance
  
  # Create NeuroVec with 4D space (spatial dims + time)
  bspace <- NeuroSpace(c(D, nobs), c(1,1,1))
  bvec <- NeuroVec(mat, bspace)
  
  # Create mask
  mask <- as.logical(NeuroVol(array(rep(1, prod(D)), D), NeuroSpace(D, c(1,1,1))))
  
  # Create labels and design
  Y <- factor(rep(c("A", "B"), length.out = nobs))
  blockVar <- rep(1:2, length.out = nobs)
  des <- mvpa_design(data.frame(Y = Y), block_var = blockVar, y_train = ~ Y)
  
  # Create mvpa_dataset
  dset <- mvpa_dataset(bvec, mask = mask)
  
  # Create a simple crossval spec
  cv_spec <- blocked_cross_validation(blockVar)
  
  # Create model spec
  mod_spec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = dset,
    design = des,
    crossval = cv_spec
  )
  
  # Calculate center voxel's linear index
  # In a 5x5x5 volume, position (3,3,3) = 3 + (3-1)*5 + (3-1)*25 = 3 + 10 + 50 = 63
  center_idx <- 63
  
  # Create a searchlight centered on the zero-variance voxel
  # Include the center and some neighbors
  vox_list <- list(
    c(center_idx - 1, center_idx, center_idx + 1)
  )
  
  # Run searchlight - should not crash when center has zero variance
  expect_no_error({
    result <- mvpa_iterate(
      mod_spec = mod_spec,
      vox_list = vox_list,
      ids = center_idx,  
      analysis_type = "searchlight",
      verbose = FALSE
    )
  })
  
  # Check that result exists and has expected structure
  expect_true(!is.null(result))
  expect_true(nrow(result) >= 1)
  
  # With the fix, the center is preserved even with zero variance
  # The analysis should complete (possibly with a warning but no crash)
})

test_that("searchlight skips ROI when center voxel has NA values", {
  library(neuroim2)
  
  D <- c(5, 5, 5)
  nobs <- 10
  
  # Create data where center voxel has NA values
  mat <- array(rnorm(prod(D) * nobs), c(D, nobs))
  center_pos <- c(3, 3, 3)
  mat[center_pos[1], center_pos[2], center_pos[3], ] <- NA  # All NA
  
  bspace <- NeuroSpace(c(D, nobs), c(1,1,1))
  bvec <- NeuroVec(mat, bspace)
  mask <- as.logical(NeuroVol(array(rep(1, prod(D)), D), NeuroSpace(D, c(1,1,1))))
  
  Y <- factor(rep(c("A", "B"), length.out = nobs))
  blockVar <- rep(1:2, length.out = nobs)
  des <- mvpa_design(data.frame(Y = Y), block_var = blockVar, y_train = ~ Y)
  dset <- mvpa_dataset(bvec, mask = mask)
  
  cv_spec <- blocked_cross_validation(blockVar)
  mod_spec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = dset,
    design = des,
    crossval = cv_spec
  )
  
  center_idx <- 63
  vox_list <- list(
    c(center_idx - 1, center_idx, center_idx + 1)
  )
  
  # Should handle gracefully without crashing
  expect_no_error({
    result <- mvpa_iterate(
      mod_spec = mod_spec,
      vox_list = vox_list,
      ids = center_idx,
      analysis_type = "searchlight",
      verbose = FALSE
    )
  })
  
  expect_true(!is.null(result))
  # The result should exist - even with NA center, it should skip gracefully
})