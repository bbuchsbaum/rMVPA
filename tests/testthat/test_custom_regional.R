library(testthat)
library(rMVPA)
library(dplyr)
library(tibble)
library(neuroim2)

# Helper function to suppress package version warnings
suppress_package_warnings <- function(expr) {
  suppressWarnings(expr, classes = "packageStartupMessage")
  # Alternative: suppress all warnings related to package versions
  # Use this if the above doesn't work:
  # withCallingHandlers(expr, 
  #   warning = function(w) {
  #     if (grepl("was built under R version", w$message)) {
  #       invokeRestart("muffleWarning")
  #     }
  #   })
}

context("run_custom_regional")

# --- Setup --- 

# Volumetric Dataset and Mask
dset_info_vol <- gen_sample_dataset(D = c(6,6,6), nobs = 40, nlevels = 2)
dataset_vol <- dset_info_vol$dataset
design_vol <- dset_info_vol$design

mask_arr <- array(0, dim(dataset_vol$mask))
mask_arr[1:3, 1:3, 1:3] <- 1 # ROI 1
mask_arr[4:6, 1:3, 1:3] <- 2 # ROI 2
mask_arr[1:3, 4:6, 4:6] <- 3 # ROI 3
mask_arr[4:6, 4:6, 4:6] <- 4 # ROI 4 (small, 1 voxel after filtering)
mask_arr[4:6, 4:6, 1:3] <- 5 # ROI 5 (will cause error in error_func)
mask_arr[1:3, 1:3, 4:6] <- 0 # Empty region
mask_arr[1,1,1] <- 0 # Ensure ROI 1 doesn't occupy the full corner
region_mask_vol <- NeuroVol(mask_arr, space(dataset_vol$mask))

# Define Custom Functions

# Returns a named list
stats_func <- function(roi_data, roi_info) {
  list(
    roi_id = roi_info$id,
    mean_val = mean(roi_data, na.rm = TRUE),
    sd_val = sd(roi_data, na.rm = TRUE),
    n_vox = ncol(roi_data)
  )
}

# Returns a single-row tibble
dataframe_func <- function(roi_data, roi_info) {
  tibble::tibble(
    avg = mean(roi_data, na.rm = TRUE),
    std = sd(roi_data, na.rm = TRUE),
    voxels = ncol(roi_data)
  )
}

# Intentionally causes an error for ROI 5
error_func <- function(roi_data, roi_info) {
  if (roi_info$id == 5) {
    stop("Intentional error in ROI 5!")
  }
  list(
    mean_val = mean(roi_data, na.rm = TRUE),
    n_vox = ncol(roi_data)
  )
}

# Returns an invalid (unnamed) list
unnamed_func <- function(roi_data, roi_info) {
  list(mean(roi_data, na.rm = TRUE), sd(roi_data, na.rm = TRUE))
}

# Returns a list with a non-scalar element
non_scalar_func <- function(roi_data, roi_info) {
  list(
    mean_val = mean(roi_data, na.rm = TRUE),
    all_vals = roi_data[,1] # Not a scalar
  )
}

# --- Basic Functionality Tests (Volumetric) --- 

test_that("run_custom_regional works with list-returning function", {
  results <- withCallingHandlers(
    run_custom_regional(dataset_vol, region_mask_vol, stats_func),
    warning = function(w) {
      if (grepl("was built under R version", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  
  expect_s3_class(results, "tbl_df")
  expect_equal(nrow(results), 5) # ROIs 1, 2, 3, 4, 5
  expect_named(results, c("id", "roi_id", "mean_val", "sd_val", "n_vox", "error", "error_message"), ignore.order = TRUE)
  expect_true(all(!results$error))
  expect_equal(results$id, results$roi_id)
  expect_true(all(results$n_vox > 0))
})

test_that("run_custom_regional works with dataframe-returning function", {
  results <- withCallingHandlers(
    run_custom_regional(dataset_vol, region_mask_vol, dataframe_func),
    warning = function(w) {
      if (grepl("was built under R version", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  
  expect_s3_class(results, "tbl_df")
  expect_equal(nrow(results), 5)
  expect_named(results, c("id", "avg", "std", "voxels", "error", "error_message"), ignore.order = TRUE)
  expect_true(all(!results$error))
  expect_true(all(results$voxels > 0))
})

# --- Error Handling within Custom Function --- 

test_that("run_custom_regional handles errors in custom_func correctly", {
  # Suppress warnings about the error itself during the run
  suppressWarnings({
      results <- withCallingHandlers(
        run_custom_regional(dataset_vol, region_mask_vol, error_func),
        warning = function(w) {
          if (grepl("was built under R version", w$message)) {
            invokeRestart("muffleWarning")
          }
        }
      )
  })
  
  expect_s3_class(results, "tbl_df")
  expect_equal(nrow(results), 5)
  expect_named(results, c("id", "mean_val", "n_vox", "error", "error_message"), ignore.order = TRUE)
  
  # Check ROI 5 (errored)
  roi5_row <- results[results$id == 5, ]
  expect_true(roi5_row$error)
  expect_match(roi5_row$error_message, "Intentional error in ROI 5!")
  expect_true(is.na(roi5_row$mean_val))
  expect_true(is.na(roi5_row$n_vox))
  
  # Check other ROIs (not errored)
  other_rows <- results[results$id != 5, ]
  expect_true(all(!other_rows$error))
  expect_true(all(!is.na(other_rows$mean_val)))
  expect_true(all(!is.na(other_rows$n_vox)))
})

# --- Parallel Execution Test --- 

test_that("run_custom_regional runs in parallel without error", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply") # Needed for auto plan setting
  
  # Run sequentially first
  results_seq <- withCallingHandlers(
    run_custom_regional(dataset_vol, region_mask_vol, stats_func, .cores = 1),
    warning = function(w) {
      if (grepl("was built under R version", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  
  # Run in parallel 
  # No need to manually set plan if future.apply is installed
  # suppressMessages to hide the plan setting message
  suppressMessages({
      results_par <- withCallingHandlers(
        run_custom_regional(dataset_vol, region_mask_vol, stats_func, .cores = 2),
        warning = function(w) {
          if (grepl("was built under R version", w$message)) {
            invokeRestart("muffleWarning")
          }
        }
      )
  })
 
  # Reset plan just in case
  future::plan(future::sequential)
  
  # Basic structural comparison (exact numeric values might differ slightly)
  expect_equal(nrow(results_par), nrow(results_seq))
  expect_equal(names(results_par), names(results_seq))
  expect_equal(results_par$id, results_seq$id)
  expect_equal(results_par$error, results_seq$error)
  expect_true(all(!results_par$error)) # Ensure no errors occurred during parallel run
})

# --- Input Validation and Edge Cases --- 

test_that("run_custom_regional input validation works", {
  plain_array <- array(0, c(2,2,2))
  # Invalid dataset
  expect_error(run_custom_regional(as.matrix(dataset_vol$train_data), region_mask_vol, stats_func),
               "`dataset` must be an 'mvpa_dataset' or 'mvpa_surface_dataset' object.")
  # Invalid mask (plain array)
  expect_error(run_custom_regional(dataset_vol, plain_array, stats_func),
               "`region_mask` must be a 'NeuroVol' or 'NeuroSurface' object.")
  # Invalid custom_func
  expect_error(run_custom_regional(dataset_vol, region_mask_vol, "not_a_function"),
               "`custom_func` must be a function.")
   # Invalid cores
  expect_error(run_custom_regional(dataset_vol, region_mask_vol, stats_func, .cores = 0),
               "`.cores` must be a positive integer.")
  expect_error(run_custom_regional(dataset_vol, region_mask_vol, stats_func, .cores = 1.5),
               "`.cores` must be a positive integer.")
})

test_that("run_custom_regional handles custom_func returning invalid structures", {
  # Unnamed list - should cause error during internal processing
  # Should result in error = TRUE in the output table
  results_unnamed <- withCallingHandlers(
    run_custom_regional(dataset_vol, region_mask_vol, unnamed_func),
    warning = function(w) {
      if (grepl("was built under R version", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_s3_class(results_unnamed, "tbl_df")
  expect_true(all(results_unnamed$error))
  expect_match(results_unnamed$error_message[1], "Error in custom_func: The list or data frame returned by custom_func must have names")
  
  # Non-scalar list - should cause error during internal processing
  # Should result in error = TRUE in the output table
  results_nonscalar <- withCallingHandlers(
    run_custom_regional(dataset_vol, region_mask_vol, non_scalar_func),
    warning = function(w) {
      if (grepl("was built under R version", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  expect_s3_class(results_nonscalar, "tbl_df")
  expect_true(all(results_nonscalar$error))
  expect_match(results_nonscalar$error_message[1], "Error in custom_func: custom_func must return a named list of scalars")
})

# --- Verbose Option --- 

test_that("run_custom_regional runs with .verbose = TRUE", {
  # Just check that it runs without error, capturing output is complex
  # Remove subsetting of mask
  expect_silent({ 
      capture.output(
        suppressWarnings(
          run_custom_regional(dataset_vol, region_mask_vol, stats_func, .verbose = TRUE)
        )
      )
  })
})

# --- Cleanup --- 
# Ensure future plan is reset
future::plan(future::sequential) 
