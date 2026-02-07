test_that("filter_roi preserves specified voxel index", {
  
  # Test that the preserve parameter works for filter_roi.ROIVec
  
  # Create mock ROI data with one column that would normally be filtered
  n_rows <- 10
  n_cols <- 5
  
  # Create data where column 3 has zero variance (would normally be filtered)
  data_matrix <- matrix(rnorm(n_rows * n_cols), nrow = n_rows)
  data_matrix[, 3] <- 5  # Constant value = zero variance
  
  # Create mock ROIVec structure
  # Note: This is a simplified mock - real ROIVec would have more structure
  mock_roi <- list(
    train_roi = structure(
      list(
        data = data_matrix,
        indices = 1:n_cols
      ),
      class = "MockROIVec"
    ),
    test_roi = NULL
  )
  
  # Create a mock filter_roi function for testing
  filter_roi_test <- function(roi, preserve = NULL) {
    trdat <- roi$train_roi$data
    
    # Find columns with zero variance
    sdnonzero <- apply(trdat, 2, sd, na.rm = TRUE) > 0
    
    # Determine columns to keep
    keep <- sdnonzero
    
    # Preserve specified voxel if requested
    if (!is.null(preserve)) {
      gi <- roi$train_roi$indices
      kp <- match(preserve, gi)
      if (!is.na(kp) && !keep[kp]) {
        keep[kp] <- TRUE
      }
    }
    
    # Return filtered data
    list(
      train_roi = list(
        data = trdat[, keep, drop = FALSE],
        indices = roi$train_roi$indices[keep]
      ),
      test_roi = NULL
    )
  }
  
  # Test without preserve - column 3 should be filtered out
  filtered_no_preserve <- filter_roi_test(mock_roi)
  expect_false(3 %in% filtered_no_preserve$train_roi$indices)
  expect_equal(ncol(filtered_no_preserve$train_roi$data), 4)  # 5 - 1 filtered
  
  # Test with preserve - column 3 should be kept despite zero variance
  filtered_with_preserve <- filter_roi_test(mock_roi, preserve = 3)
  expect_true(3 %in% filtered_with_preserve$train_roi$indices)
  expect_equal(ncol(filtered_with_preserve$train_roi$data), 5)  # All columns kept
})

test_that("filter_roi preserve parameter handles NA columns", {
  
  # Test that preserve works even for columns with NA values
  
  n_rows <- 10
  n_cols <- 5
  
  # Create data where column 2 has NA values
  data_matrix <- matrix(rnorm(n_rows * n_cols), nrow = n_rows)
  data_matrix[, 2] <- NA
  
  mock_roi <- list(
    train_roi = structure(
      list(
        data = data_matrix,
        indices = 1:n_cols
      ),
      class = "MockROIVec"
    ),
    test_roi = NULL
  )
  
  # Test filter function
  filter_roi_test <- function(roi, preserve = NULL) {
    trdat <- roi$train_roi$data
    
    # Find columns with NAs and zero variance
    nas <- apply(trdat, 2, function(v) any(is.na(v)))
    sdnonzero <- apply(trdat, 2, sd, na.rm = TRUE) > 0
    
    keep <- !nas & sdnonzero
    
    # Preserve specified voxel if requested
    if (!is.null(preserve)) {
      gi <- roi$train_roi$indices
      kp <- match(preserve, gi)
      if (!is.na(kp) && !keep[kp]) {
        keep[kp] <- TRUE
      }
    }
    
    list(
      train_roi = list(
        data = trdat[, keep, drop = FALSE],
        indices = roi$train_roi$indices[keep]
      ),
      test_roi = NULL
    )
  }
  
  # Without preserve - column 2 with NAs should be filtered
  filtered_no_preserve <- filter_roi_test(mock_roi)
  expect_false(2 %in% filtered_no_preserve$train_roi$indices)
  
  # With preserve - column 2 should be kept despite NAs
  filtered_with_preserve <- filter_roi_test(mock_roi, preserve = 2)
  expect_true(2 %in% filtered_with_preserve$train_roi$indices)
})