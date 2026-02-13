library(testthat)
library(rMVPA)

# Skip entire file if shard is not installed
skip_if_not_installed("shard")

test_that("use_shard tags model spec with shard_model_spec class", {
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  mspec_shard <- use_shard(mspec)
  expect_true(inherits(mspec_shard, "shard_model_spec"))
  expect_false(is.null(mspec_shard$shard_data))
  expect_false(is.null(mspec_shard$shard_data$shared_train))
  expect_false(is.null(mspec_shard$shard_data$idx_to_col))
  expect_equal(mspec_shard$shard_data$roi_type, "volumetric")

  shard_cleanup(mspec_shard$shard_data)
})

test_that("shard_extract_roi produces valid ROIVec", {
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec_shard <- use_shard(mspec)

  # Get some voxel indices from the mask
  mask_idx <- which(as.logical(mspec_shard$dataset$mask))
  vox <- mask_idx[1:min(10, length(mask_idx))]

  roi <- rMVPA:::shard_extract_roi(vox, mspec_shard$shard_data,
                                    min_voxels = 1)
  expect_false(is.null(roi))
  expect_true(inherits(roi$train_roi, "ROIVec"))
  expect_equal(nrow(neuroim2::values(roi$train_roi)), 20)

  shard_cleanup(mspec_shard$shard_data)
})

test_that("shard searchlight produces same results as default", {
  skip_on_cran()

  ds <- gen_sample_dataset(c(5,5,5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  # Get a small set of searchlight spheres
  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  # Run default path
  set.seed(42)
  res_default <- mvpa_iterate(mspec, vox_iter,
                               ids = seq_along(vox_iter))

  # Run shard path
  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_shard <- mvpa_iterate(mspec_shard, vox_iter,
                             ids = seq_along(vox_iter))

  # Both should produce results with same number of rows
  expect_equal(nrow(res_default), nrow(res_shard))

  # Same IDs
  expect_equal(sort(res_default$id), sort(res_shard$id))

  # Same error status
  expect_equal(res_default$error, res_shard$error)

  # Deeper parity: compare performance values (audit #3)
  if (!all(res_default$error)) {
    perf_default <- unlist(lapply(res_default$performance, function(p) {
      if (is.null(p)) NA_real_ else p[[1]]
    }))
    perf_shard <- unlist(lapply(res_shard$performance, function(p) {
      if (is.null(p)) NA_real_ else p[[1]]
    }))
    expect_equal(perf_default, perf_shard)
  }
})

# ---- Audit #1: shard_cleanup uses close(), not unshare() ----

test_that("shard_cleanup runs without error", {
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec_shard <- use_shard(mspec)

  # Cleanup should not error

  expect_no_error(shard_cleanup(mspec_shard$shard_data))
  # Repeated cleanup on the same handles should also be safe
  expect_no_error(shard_cleanup(mspec_shard$shard_data))

  # Cleanup on NULL is safe
  expect_no_error(shard_cleanup(NULL))
})

# ---- Audit #2: multibasis datasets are rejected early ----

test_that("use_shard rejects multibasis datasets with informative error", {
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks = 2, nlevels = 2)
  mb <- mvpa_multibasis_dataset(
    train_data = list(ds$dataset$train_data, ds$dataset$train_data),
    mask = ds$dataset$mask
  )
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, mb, ds$design,
                      "classification", crossval = cval)

  expect_error(use_shard(mspec), "multibasis.*not yet supported")
})

# ---- Audit #4: repeated use_shard() cleans up old handles ----

test_that("repeated use_shard() cleans up previous handles", {
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  mspec1 <- use_shard(mspec)
  old_data <- mspec1$shard_data

  # Second call should not error (cleans up old handles)
  mspec2 <- use_shard(mspec1)

  expect_true(inherits(mspec2, "shard_model_spec"))
  expect_false(is.null(mspec2$shard_data))
  # Class should not be duplicated
  expect_equal(sum(class(mspec2) == "shard_model_spec"), 1L)

  shard_cleanup(mspec2$shard_data)
})

# ---- Audit #3: cross-process ALTREP serialization via multisession ----

test_that("shard backend works across PSOCK multisession workers", {
  skip_on_cran()

  ds <- gen_sample_dataset(c(5,5,5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  # Force PSOCK multisession (2 workers) to test ALTREP serialization
  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)

  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_shard <- mvpa_iterate(mspec_shard, vox_iter,
                             ids = seq_along(vox_iter))

  expect_equal(nrow(res_shard), length(vox_iter))
  expect_equal(sort(res_shard$id), sort(seq_along(vox_iter)))
  # No errors from deserialization failures
  expect_true(all(!res_shard$error))
})

test_that("shard surface searchlight parity with default backend", {
  skip_on_cran()
  skip_if_not_installed("neurosurf")

  ds <- gen_sample_dataset(
    D = 100,
    nobs = 20,
    blocks = 2,
    nlevels = 2,
    data_mode = "surface"
  )
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 8)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  set.seed(42)
  res_default <- mvpa_iterate(
    mspec, vox_iter,
    ids = seq_along(vox_iter),
    analysis_type = "searchlight"
  )

  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_shard <- mvpa_iterate(
    mspec_shard, vox_iter,
    ids = seq_along(vox_iter),
    analysis_type = "searchlight"
  )

  expect_equal(nrow(res_default), nrow(res_shard))
  expect_equal(sort(res_default$id), sort(res_shard$id))
  expect_equal(res_default$error, res_shard$error)
})

test_that("shard clustered regional parity with default backend", {
  skip_on_cran()

  ds <- gen_clustered_sample_dataset(
    D = c(8, 8, 8),
    nobs = 20,
    K = 8,
    nlevels = 2,
    blocks = 2
  )
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, k = 4)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  set.seed(42)
  res_default <- mvpa_iterate(
    mspec, vox_iter,
    ids = seq_along(vox_iter),
    analysis_type = "regional"
  )

  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_shard <- mvpa_iterate(
    mspec_shard, vox_iter,
    ids = seq_along(vox_iter),
    analysis_type = "regional"
  )

  expect_equal(nrow(res_default), nrow(res_shard))
  expect_equal(sort(res_default$id), sort(res_shard$id))
  expect_equal(res_default$error, res_shard$error)
})

test_that("shard parity with external test set", {
  skip_on_cran()

  ds <- gen_sample_dataset(
    c(5, 5, 5), 20,
    blocks = 2,
    nlevels = 2,
    external_test = TRUE
  )
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  set.seed(42)
  res_default <- mvpa_iterate(mspec, vox_iter, ids = seq_along(vox_iter))

  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_shard <- mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))

  expect_equal(nrow(res_default), nrow(res_shard))
  expect_equal(sort(res_default$id), sort(res_shard$id))
  expect_equal(res_default$error, res_shard$error)
})

test_that("shard_extract_roi drops invalid indices cleanly", {
  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec_shard <- use_shard(mspec)

  valid_vox <- which(as.logical(mspec_shard$dataset$mask))[1]
  bad_vox <- c(
    valid_vox,
    0L,
    NA_integer_,
    max(which(as.logical(mspec_shard$dataset$mask))) + 1000L
  )

  roi <- rMVPA:::shard_extract_roi(bad_vox, mspec_shard$shard_data, min_voxels = 1)

  expect_false(is.null(roi))
  expect_equal(length(neuroim2::indices(roi$train_roi)), 1L)
  expect_equal(as.integer(neuroim2::indices(roi$train_roi)), as.integer(valid_vox))

  shard_cleanup(mspec_shard$shard_data)
})

test_that("post-run shard cleanup is idempotent", {
  skip_on_cran()

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec_shard <- use_shard(mspec)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(3, length(vox_iter))]

  expect_no_error(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )
  # mvpa_iterate() already triggers shard_cleanup() via on.exit
  expect_no_error(shard_cleanup(mspec_shard$shard_data))
  expect_no_error(shard_cleanup(mspec_shard$shard_data))
})

test_that("run_searchlight backend='shard' works without explicit use_shard", {
  skip_on_cran()

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  expect_no_error({
    res <- run_searchlight(
      mspec,
      radius = 3,
      method = "standard",
      backend = "shard",
      batch_size = 5
    )
    expect_true(inherits(res, "searchlight_result"))
    expect_true(length(res$results) > 0)
  })
})

test_that("run_regional backend='shard' works without explicit use_shard", {
  skip_on_cran()

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  region_mask <- neuroim2::NeuroVol(
    sample(1:4, length(ds$dataset$mask), replace = TRUE),
    neuroim2::space(ds$dataset$mask)
  )

  expect_no_error({
    res <- run_regional(mspec, region_mask, backend = "shard")
    expect_true(inherits(res, "regional_mvpa_result"))
    expect_true(is.data.frame(res$performance_table))
  })
})

test_that("backend='auto' runs and remains backward compatible", {
  skip_on_cran()

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  expect_no_error({
    res <- run_searchlight(
      mspec,
      radius = 3,
      method = "standard",
      backend = "auto",
      batch_size = 5
    )
    expect_true(inherits(res, "searchlight_result"))
  })
})

test_that("shard workers do not inherit large caller globals", {
  skip_on_cran()
  skip_if_not_installed("future")

  # Keep globals budget small enough that accidental large exports fail fast.
  withr::local_options(list(future.globals.maxSize = 8 * 1024^2))
  huge_caller_object <- rnorm(2 * 1024^2) # ~16 MB

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec <- use_shard(mspec)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(3, length(vox_iter))]

  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)

  expect_no_error(
    mvpa_iterate(mspec, vox_iter, ids = seq_along(vox_iter))
  )

  expect_gt(length(huge_caller_object), 0L)
})

test_that("shard backend uses carrier::crate when carrier is installed", {
  skip_on_cran()
  skip_if_not_installed("carrier")

  withr::local_options(list(rMVPA.test.used_carrier_crate = FALSE))
  trace("crate",
        where = asNamespace("carrier"),
        tracer = quote(options(rMVPA.test.used_carrier_crate = TRUE)),
        print = FALSE)
  on.exit(untrace("crate", where = asNamespace("carrier")), add = TRUE)

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec <- use_shard(mspec)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(2, length(vox_iter))]

  expect_no_error(
    mvpa_iterate(mspec, vox_iter, ids = seq_along(vox_iter))
  )
  expect_true(isTRUE(getOption("rMVPA.test.used_carrier_crate")))
})
