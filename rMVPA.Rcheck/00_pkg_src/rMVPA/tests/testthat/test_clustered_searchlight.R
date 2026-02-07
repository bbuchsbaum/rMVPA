testthat::skip_if_not_installed("neuroim2")
library(neuroim2)

test_that("mvpa_clustered_dataset constructor works", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 5)
  dset <- ds$dataset

  expect_s3_class(dset, "mvpa_clustered_dataset")
  expect_s3_class(dset, "mvpa_dataset")
  expect_equal(nobs(dset), 20L)
  expect_false(has_test_set(dset))
  expect_true(inherits(dset$cvol, "ClusteredNeuroVol"))
  expect_true(inherits(dset$region_mask, "NeuroVol"))
})

test_that("mvpa_clustered_dataset with test set", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 5,
                                      external_test = TRUE)
  dset <- ds$dataset

  expect_true(has_test_set(dset))
  expect_true(!is.null(dset$test_data))
})

test_that("get_center_ids returns cluster IDs", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 5)
  cids <- get_center_ids(ds$dataset)
  K <- num_clusters(ds$dataset$cvol)

  expect_equal(cids, seq_len(K))
  expect_true(K >= 2)
})

test_that("searchlight_scope returns 'regional' for clustered dataset", {
  ds <- gen_clustered_sample_dataset()
  expect_equal(rMVPA:::searchlight_scope(ds$dataset), "regional")
})

test_that("get_searchlight returns correct specs with k-nearest", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 5)
  sl <- get_searchlight(ds$dataset, "standard", k = 3)
  K <- num_clusters(ds$dataset$cvol)

  expect_equal(length(sl), K)
  # Each spec should have at most k+1 neighbors (including self)
  for (spec in sl) {
    expect_s3_class(spec, "clustered_roi_spec")
    expect_true(length(spec$neighbors) >= 1)
    expect_true(length(spec$neighbors) <= 4)  # k+1 = 4
    # Seed should be in neighbors
    expect_true(spec$seed %in% spec$neighbors)
  }
})

test_that("get_searchlight returns correct specs with radius", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 5)
  # Use a large radius to include all clusters
  sl <- get_searchlight(ds$dataset, "standard", radius = 100, k = NULL)
  K <- num_clusters(ds$dataset$cvol)

  expect_equal(length(sl), K)
  # With radius=100 on a 10x10x10 volume, all clusters should be neighbors
  for (spec in sl) {
    expect_equal(length(spec$neighbors), K)
  }
})

test_that("get_samples and as_roi work for clustered datasets", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 5)
  sl <- get_searchlight(ds$dataset, "standard", k = 3)

  samples <- get_samples(ds$dataset, sl)
  expect_equal(nrow(samples), length(sl))
  expect_true("sample" %in% names(samples))

  # Test as_roi for first sample
  sam <- samples$sample[[1]]
  roi <- as_roi(sam, ds$dataset)

  expect_true(is.list(roi))
  expect_true("train_roi" %in% names(roi))
  expect_true(inherits(roi$train_roi, "ROIVec"))
  expect_null(roi$test_roi)
  # Rows should be timepoints, columns should be neighbors
  vals <- neuroim2::values(roi$train_roi)
  expect_equal(nrow(vals), 20)
  expect_equal(ncol(vals), length(sam$neighbors))
})

test_that("build_output_map broadcasts to all voxels in cluster", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 20, K = 5)
  K <- num_clusters(ds$dataset$cvol)
  ids <- seq_len(K)
  vals <- seq_len(K) * 10.0

  out_map <- rMVPA:::build_output_map(ds$dataset, vals, ids)
  expect_true(inherits(out_map, "NeuroVol"))

  # Check that each cluster's voxels got the right value
  rmask <- ds$dataset$region_mask
  for (k in seq_len(K)) {
    vox_in_cluster <- which(as.array(rmask) == k)
    if (length(vox_in_cluster) > 0) {
      expect_true(all(out_map[vox_in_cluster] == k * 10.0))
    }
  }
})

test_that("full clustered searchlight pipeline runs with sda_notune", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 5,
                                      nlevels = 2, blocks = 3)
  cval <- blocked_cross_validation(ds$design$block_var)
  model <- load_model("sda_notune")

  mspec <- mvpa_model(
    model      = model,
    dataset    = ds$dataset,
    design     = ds$design,
    model_type = "classification",
    crossval   = cval
  )

  results <- run_searchlight(mspec, radius = NULL, method = "standard", k = 3)

  expect_true(inherits(results, "searchlight_result"))
  expect_true(length(results$results) > 0)
  expect_true(results$active_voxels > 0)
})

test_that("full clustered searchlight pipeline with radius", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 5,
                                      nlevels = 2, blocks = 3)
  cval <- blocked_cross_validation(ds$design$block_var)
  model <- load_model("sda_notune")

  mspec <- mvpa_model(
    model      = model,
    dataset    = ds$dataset,
    design     = ds$design,
    model_type = "classification",
    crossval   = cval
  )

  results <- run_searchlight(mspec, radius = 20, method = "standard", k = NULL)

  expect_true(inherits(results, "searchlight_result"))
  expect_true(length(results$results) > 0)
})

test_that("clustered searchlight with external test set", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 5,
                                      nlevels = 2, blocks = 3,
                                      external_test = TRUE)
  model <- load_model("sda_notune")

  mspec <- mvpa_model(
    model      = model,
    dataset    = ds$dataset,
    design     = ds$design,
    model_type = "classification"
  )

  results <- run_searchlight(mspec, radius = NULL, method = "standard", k = 3)

  expect_true(inherits(results, "searchlight_result"))
  expect_true(length(results$results) > 0)
})

test_that("clustered searchlight with small k works", {
  ds <- gen_clustered_sample_dataset(D = c(10, 10, 10), nobs = 30, K = 5,
                                      nlevels = 2, blocks = 3)
  cval <- blocked_cross_validation(ds$design$block_var)
  model <- load_model("sda_notune")

  mspec <- mvpa_model(
    model      = model,
    dataset    = ds$dataset,
    design     = ds$design,
    model_type = "classification",
    crossval   = cval
  )

  # k=2 means each parcel uses itself + 2 nearest
  results <- run_searchlight(mspec, radius = NULL, method = "standard", k = 2)

  expect_true(inherits(results, "searchlight_result"))
  expect_true(length(results$results) > 0)
})

test_that("print method works for mvpa_clustered_dataset", {
  ds <- gen_clustered_sample_dataset()
  expect_output(print(ds$dataset), "Clustered MVPA Dataset")
})

test_that("length.clustered_roi_spec returns neighbor count", {
  spec <- rMVPA:::clustered_roi_spec(seed = 1L, neighbors = c(1L, 2L, 3L))
  expect_equal(length(spec), 3)
})
