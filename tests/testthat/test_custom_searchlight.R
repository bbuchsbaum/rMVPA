library(testthat)
library(rMVPA)
library(neuroim2)

context("run_custom_searchlight")

# --- Setup ---

# Generate a sample volumetric dataset
dset_info_vol <- gen_sample_dataset(D = c(5, 5, 5), nobs = 20, nlevels = 2)
dataset_vol <- dset_info_vol$dataset

# Define a simple custom function for the searchlight
# Assume it receives data (samples x voxels_in_sphere) and info
# It should return a single named value or a list with one named value
mean_signal_sl <- function(sl_data, sl_info) {
  # sl_data: matrix of samples x voxels within the sphere
  # sl_info: list containing info like center voxel index, coords etc.
  #         (Exact structure depends on the final implementation)
  mean_val <- mean(sl_data, na.rm = TRUE)
  # Return a named list with one scalar value
  list(mean_signal = mean_val)
}

# --- Basic Functionality Test ---

test_that("run_custom_searchlight (standard) runs without error and returns correct structure", {
  # Run standard searchlight
  searchlight_results <- run_custom_searchlight(
    dataset = dataset_vol,
    custom_func = mean_signal_sl,
    radius = 5, # Use a slightly smaller radius for faster testing
    method = "standard",
    .cores = 1, # Keep it simple first
    .verbose = FALSE
  )

  # Check main object class
  expect_s3_class(searchlight_results, "searchlight_result")
  expect_true(is.list(searchlight_results))
  expect_named(searchlight_results, c("results", "n_voxels", "active_voxels", "metrics"))

  # Check metrics list
  expect_equal(searchlight_results$metrics, c("mean_signal"))
  expect_true(is.list(searchlight_results$results))
  expect_named(searchlight_results$results, c("mean_signal"))

  # Check the performance object for the metric
  perf_obj <- searchlight_results$results$mean_signal
  expect_s3_class(perf_obj, "searchlight_performance")
  expect_named(perf_obj, c("data", "metric_name", "n_nonzero", "summary_stats", "indices"))

  # Check the actual data map (NeuroVol)
  map_vol <- perf_obj$data
  expect_s3_class(map_vol, "NeuroVol")
  expect_equal(dim(map_vol), dim(dataset_vol$mask))
  expect_equal(space(map_vol), space(dataset_vol$mask))
  expect_true(is.numeric(values(map_vol)))
  
  # Check that some valid (non-NA) results were computed in the active mask areas
  active_indices <- which(as.logical(dataset_vol$mask))
  expect_false(all(is.na(values(map_vol)[active_indices]))) 
  
  # Check summary stats are populated
   expect_true(is.list(perf_obj$summary_stats))
   expect_named(perf_obj$summary_stats, c("mean", "sd", "min", "max"))
   expect_true(all(sapply(perf_obj$summary_stats, is.numeric))) 
   
   # Check indices (should be center voxels for standard)
    expect_true(is.numeric(perf_obj$indices))
    expect_equal(sort(perf_obj$indices), sort(active_indices)) # Standard covers all active centers
    
})


test_that("run_custom_searchlight (randomized) runs without error", {
  skip_on_cran() # May take a bit longer
  
  # Run randomized searchlight
  searchlight_results_rand <- run_custom_searchlight(
    dataset = dataset_vol,
    custom_func = mean_signal_sl,
    radius = 5,
    method = "randomized",
    niter = 10, # Fewer iterations for testing
    .cores = 1,
    .verbose = FALSE
  )

  # Basic structure checks (similar to standard)
  expect_s3_class(searchlight_results_rand, "searchlight_result")
  expect_named(searchlight_results_rand, c("results", "n_voxels", "active_voxels", "metrics"))
  expect_equal(searchlight_results_rand$metrics, c("mean_signal"))
  expect_s3_class(searchlight_results_rand$results$mean_signal, "searchlight_performance")
  map_vol_rand <- searchlight_results_rand$results$mean_signal$data
  expect_s3_class(map_vol_rand, "NeuroVol")
  expect_equal(dim(map_vol_rand), dim(dataset_vol$mask))
  
  # Check that some results exist (might not cover all voxels unlike standard)
  active_indices <- which(as.logical(dataset_vol$mask))
  expect_false(all(is.na(values(map_vol_rand)[active_indices]))) 
  
   # Indices should be NULL for randomized combined results
    expect_null(searchlight_results_rand$results$mean_signal$indices)
})


test_that("run_custom_searchlight handles errors in custom_func", {
  # Define a function that errors if fewer than 3 voxels in sphere
  error_sl_func <- function(sl_data, sl_info) {
    if (ncol(sl_data) < 3) {
      stop("Test Error: Too few voxels!")
    }
    list(mean_signal = mean(sl_data, na.rm = TRUE))
  }

  # Run with a small radius likely to trigger the error
  # Suppress warnings expected from the error handling during run
  suppressWarnings({
      searchlight_results_err <- run_custom_searchlight(
          dataset = dataset_vol,
          custom_func = error_sl_func,
          radius = 1, # Very small radius
          method = "standard",
          .cores = 1,
          .verbose = FALSE
      )
  })

  # Check structure is still valid
  expect_s3_class(searchlight_results_err, "searchlight_result")
  expect_named(searchlight_results_err$results, "mean_signal")
  map_vol_err <- searchlight_results_err$results$mean_signal$data
  expect_s3_class(map_vol_err, "NeuroVol")

  # Check that some values are NA (where the error occurred)
  # but not necessarily all, depending on sphere sizes
  active_indices <- which(as.logical(dataset_vol$mask))
  map_values_active <- values(map_vol_err)[active_indices]
  expect_true(any(is.na(map_values_active)))
  # Ensure not ALL are NA (unless the mask is tiny/radius catches none)
  if (length(active_indices) > 0) {
      expect_false(all(is.na(map_values_active))) 
  }
})

test_that("run_custom_searchlight runs in parallel (standard)", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  # Run sequentially
  results_seq <- run_custom_searchlight(
    dataset = dataset_vol,
    custom_func = mean_signal_sl,
    radius = 5,
    method = "standard",
    .cores = 1, .verbose = FALSE
  )

  # Run in parallel
  suppressMessages({
      results_par <- run_custom_searchlight(
        dataset = dataset_vol,
        custom_func = mean_signal_sl,
        radius = 5,
        method = "standard",
        .cores = 2, .verbose = FALSE
      )
  })
  
  # Reset plan
  future::plan(future::sequential)

  # Compare structure and basic properties
  expect_equal(names(results_par), names(results_seq))
  expect_equal(results_par$metrics, results_seq$metrics)
  expect_equal(dim(results_par$results$mean_signal$data), 
               dim(results_seq$results$mean_signal$data))
               
  # Compare numeric results (should be identical for standard method)
   expect_equal(values(results_par$results$mean_signal$data),
                values(results_seq$results$mean_signal$data)) 

})
