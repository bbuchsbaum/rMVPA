test_that("save_results defaults to individual files", {
  library(neuroim2)
  
  # Create mock searchlight result with multiple metrics
  dim <- c(10, 10, 10)
  space <- NeuroSpace(dim, c(1,1,1))
  
  # Create three different metric volumes
  accuracy_vol <- NeuroVol(array(runif(prod(dim), 0.5, 1.0), dim), space)
  auc_vol <- NeuroVol(array(runif(prod(dim), 0.4, 0.9), dim), space)
  sensitivity_vol <- NeuroVol(array(runif(prod(dim), 0.3, 0.8), dim), space)
  
  # Create searchlight_result object
  result <- structure(
    list(
      results = list(
        accuracy = accuracy_vol,
        auc = auc_vol,
        sensitivity = sensitivity_vol
      ),
      n_voxels = prod(dim),
      active_voxels = prod(dim),
      metrics = c("accuracy", "auc", "sensitivity")
    ),
    class = c("searchlight_result", "list")
  )
  
  # Create temp directory for output
  temp_dir <- tempfile()
  
  # Save with defaults (should create individual files)
  save_results(result, temp_dir, quiet = TRUE)
  
  # Check that individual files were created
  map_files <- list.files(file.path(temp_dir, "maps"), pattern = "\\.nii\\.gz$")
  
  # Should have 3 separate files
  expect_equal(length(map_files), 3)
  
  # Check that files are named after metrics
  expect_true("accuracy.nii.gz" %in% map_files)
  expect_true("auc.nii.gz" %in% map_files)
  expect_true("sensitivity.nii.gz" %in% map_files)
  
  # Check manifest was created (standard level is default)
  manifest_files <- list.files(temp_dir, pattern = "^manifest\\.")
  expect_true(length(manifest_files) > 0)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("save_results can still stack into 4D when requested", {
  library(neuroim2)
  
  dim <- c(10, 10, 10)
  space <- NeuroSpace(dim, c(1,1,1))
  
  accuracy_vol <- NeuroVol(array(runif(prod(dim), 0.5, 1.0), dim), space)
  auc_vol <- NeuroVol(array(runif(prod(dim), 0.4, 0.9), dim), space)
  
  result <- structure(
    list(
      results = list(
        accuracy = accuracy_vol,
        auc = auc_vol
      ),
      n_voxels = prod(dim),
      active_voxels = prod(dim),
      metrics = c("accuracy", "auc")
    ),
    class = c("searchlight_result", "list")
  )
  
  temp_dir <- tempfile()
  
  # Save with stack="vec" to force 4D output
  save_results(result, temp_dir, stack = "vec", quiet = TRUE)
  
  # Check that only one file was created
  map_files <- list.files(file.path(temp_dir, "maps"), pattern = "\\.nii\\.gz$")
  expect_equal(length(map_files), 1)
  
  # Default name should be searchlight.nii.gz
  expect_true("searchlight.nii.gz" %in% map_files)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("save_results handles mixed metric names appropriately", {
  library(neuroim2)
  
  dim <- c(10, 10, 10)
  space <- NeuroSpace(dim, c(1,1,1))
  
  # Create volumes with various naming challenges
  vol1 <- NeuroVol(array(runif(prod(dim)), dim), space)
  vol2 <- NeuroVol(array(runif(prod(dim)), dim), space)
  vol3 <- NeuroVol(array(runif(prod(dim)), dim), space)
  
  result <- structure(
    list(
      results = list(
        `Test Accuracy (%)` = vol1,  # Spaces and special chars
        `ROC-AUC` = vol2,            # Hyphen
        sensitivity_specificity = vol3  # Underscore
      ),
      n_voxels = prod(dim),
      active_voxels = prod(dim),
      metrics = c("Test Accuracy (%)", "ROC-AUC", "sensitivity_specificity")
    ),
    class = c("searchlight_result", "list")
  )
  
  temp_dir <- tempfile()
  
  save_results(result, temp_dir, quiet = TRUE)
  
  map_files <- list.files(file.path(temp_dir, "maps"), pattern = "\\.nii\\.gz$")
  
  # Should have 3 files with slugified names
  expect_equal(length(map_files), 3)
  
  # Check that problematic characters were handled
  expect_true(any(grepl("Test_Accuracy", map_files)))
  expect_true(any(grepl("ROC-AUC", map_files)))  # Hyphens are preserved by slugify
  expect_true("sensitivity_specificity.nii.gz" %in% map_files)
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})