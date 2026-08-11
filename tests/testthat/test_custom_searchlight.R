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
    radius = 2, # Maximum radius for 5x5x5 dataset
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
  expect_true(inherits(map_vol, "NeuroVol"))
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
    # With smaller radius, some edge voxels might not have valid searchlight spheres
    expect_true(length(perf_obj$indices) <= length(active_indices))
    expect_true(all(perf_obj$indices %in% active_indices))
    
})


test_that("run_custom_searchlight (randomized) runs without error", {
  # Create a fresh dataset for this test to ensure consistency
  dset_info_rand <- gen_sample_dataset(D = c(6, 6, 6), nobs = 20, nlevels = 2)
  dataset_rand <- dset_info_rand$dataset
  
  # Run randomized searchlight
  searchlight_results_rand <- run_custom_searchlight(
    dataset = dataset_rand,
    custom_func = mean_signal_sl,
    radius = 3, # Maximum radius for 6x6x6 dataset
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
  expect_true(inherits(map_vol_rand, "NeuroVol"))
  expect_equal(dim(map_vol_rand), dim(dataset_rand$mask))
  
  # Check that some results exist (might not cover all voxels unlike standard)
  active_indices <- which(as.logical(dataset_rand$mask))
  expect_false(all(is.na(values(map_vol_rand)[active_indices]))) 
  
   # Indices should be NULL for randomized combined results
    expect_null(searchlight_results_rand$results$mean_signal$indices)
})


test_that("run_custom_searchlight handles errors in custom_func", {
  # Define a function that errors based on deterministic criteria
  # Error on specific center indices to ensure a mix of success/failure
  error_sl_func <- function(sl_data, sl_info) {
    # Use modulo arithmetic to ensure some spheres fail and some succeed
    # Error on indices where (center_index %% 4) == 1 to ensure we hit some
    if ((sl_info$center_index %% 4) == 1) {
      stop("Test Error: Selected spheres fail!")
    }
    list(mean_signal = mean(sl_data, na.rm = TRUE))
  }

  # Run with standard searchlight
  # Suppress warnings expected from the error handling during run
  suppressWarnings({
      searchlight_results_err <- run_custom_searchlight(
          dataset = dataset_vol,  # Use the original dataset
          custom_func = error_sl_func,
          radius = 2, # Maximum radius for 5x5x5 dataset
          method = "standard",
          .cores = 1,
          .verbose = FALSE
      )
  })

  # Check structure is still valid
  expect_s3_class(searchlight_results_err, "searchlight_result")
  expect_named(searchlight_results_err$results, "mean_signal")
  map_vol_err <- searchlight_results_err$results$mean_signal$data
  expect_true(inherits(map_vol_err, "NeuroVol"))

  # Get values from the result map
  all_values <- values(map_vol_err)
  active_indices <- which(as.logical(dataset_vol$mask))
  active_values <- all_values[active_indices]
  
  # Check that the searchlight completed despite errors
  # The current implementation skips failed spheres rather than inserting NAs
  # So we check that we got a valid result with fewer processed spheres
  expect_true(searchlight_results_err$active_voxels < searchlight_results_err$n_voxels,
              info = "Expected fewer active voxels due to failed spheres")
  
  # The performance object should still be valid
  expect_true(is.numeric(searchlight_results_err$results$mean_signal$n_nonzero))
  expect_true(searchlight_results_err$results$mean_signal$n_nonzero > 0)
})

test_that("run_custom_searchlight runs in parallel (standard)", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  # Run sequentially
  results_seq <- run_custom_searchlight(
    dataset = dataset_vol,
    custom_func = mean_signal_sl,
    radius = 2, # Maximum radius for 5x5x5 dataset
    method = "standard",
    .cores = 1, .verbose = FALSE
  )

  # Run in parallel
  suppressMessages({
      results_par <- run_custom_searchlight(
        dataset = dataset_vol,
        custom_func = mean_signal_sl,
        radius = 2, # Maximum radius for 5x5x5 dataset
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

test_that("custom callbacks receive separate test data and arbitrary user data", {
  set.seed(42)
  dims <- c(4, 4, 4)
  n_obs <- 6L
  sp <- neuroim2::NeuroSpace(c(dims, n_obs))
  mask_sp <- neuroim2::NeuroSpace(dims)
  covariate <- c(-2, -1, 0, 1, 2, NA_real_)

  train_array <- array(rnorm(prod(dims) * n_obs), c(dims, n_obs))
  test_array <- train_array
  for (i in seq_len(n_obs)) {
    test_array[, , , i] <- test_array[, , , i] + ifelse(is.na(covariate[i]), 0, covariate[i])
  }

  dataset <- mvpa_dataset(
    neuroim2::NeuroVec(train_array, sp),
    test_data = neuroim2::NeuroVec(test_array, sp),
    mask = neuroim2::NeuroVol(array(1, dims), mask_sp)
  )

  continuous_association <- function(sl_data, sl_info) {
    expect_true(is.matrix(sl_info$test_data))
    effect <- rowMeans(sl_info$test_data - sl_data)
    ok <- is.finite(effect) & is.finite(sl_info$user_data$covariate)
    list(
      covariate_cor = stats::cor(effect[ok], sl_info$user_data$covariate[ok]),
      n_pairs = sum(ok)
    )
  }

  result <- run_custom_searchlight(
    dataset,
    continuous_association,
    radius = 1,
    user_data = list(covariate = covariate),
    batch_size = 8,
    .cores = 1,
    .verbose = FALSE
  )

  active <- which(dataset$mask > 0)
  cors <- neuroim2::values(result$results$covariate_cor$data)[active]
  counts <- neuroim2::values(result$results$n_pairs$data)[active]
  expect_equal(cors[is.finite(cors)], rep(1, sum(is.finite(cors))), tolerance = 1e-12)
  expect_equal(counts[is.finite(counts)], rep(5, sum(is.finite(counts))))
})

test_that(".cores = 1 overrides and restores an inherited parallel plan", {
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")
  skip_on_cran()

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2)
  parent_pid <- Sys.getpid()

  result <- run_custom_searchlight(
    dataset_vol,
    function(sl_data, sl_info) list(pid = Sys.getpid()),
    radius = 1,
    batch_size = 8,
    .cores = 1,
    .verbose = FALSE
  )

  active <- which(dataset_vol$mask > 0)
  observed <- unique(stats::na.omit(
    as.numeric(neuroim2::values(result$results$pid$data))[active]
  ))
  expect_equal(observed, parent_pid)
  expect_equal(future::nbrOfWorkers(), 2L)
})

test_that("save_results writes custom searchlight performance maps", {
  result <- run_custom_searchlight(
    dataset_vol,
    mean_signal_sl,
    radius = 1,
    batch_size = 8,
    .cores = 1,
    .verbose = FALSE
  )
  out <- tempfile("custom-searchlight-save-")

  paths <- save_results(result, out, stack = "none", quiet = TRUE)

  expect_true(dir.exists(file.path(out, "maps")))
  expect_true(file.exists(file.path(out, "maps", "mean_signal.nii.gz")))
  expect_false(dir.exists(file.path(out, "aux")))
  expect_true(length(paths$maps) == 1L)
})
