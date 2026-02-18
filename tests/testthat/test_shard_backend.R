library(testthat)
library(rMVPA)

# Skip entire file if shard is not installed
skip_if_not_installed("shard")

# Muffle environment-specific package version warnings emitted on some workers
# (e.g., "package 'purrr' was built under R version ..."), while preserving
# all other warnings.
muffle_worker_version_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("was built under R version", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

run_shard_with_plan <- function(mspec, vox_iter, strategy, workers = NULL) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  if (is.null(workers)) {
    future::plan(strategy)
  } else {
    future::plan(strategy, workers = workers)
  }

  mspec_shard <- use_shard(mspec)
  on.exit(shard_cleanup(mspec_shard$shard_data), add = TRUE)

  set.seed(42)
  muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )
}

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
  res_default <- muffle_worker_version_warnings(
    mvpa_iterate(mspec, vox_iter, ids = seq_along(vox_iter))
  )

  # Run shard path
  mspec_shard <- use_shard(mspec)
  set.seed(42)
  res_shard <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )

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
  skip_if_not_installed("future")

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
  res_shard <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )

  expect_equal(nrow(res_shard), length(vox_iter))
  expect_equal(sort(res_shard$id), sort(seq_along(vox_iter)))
  # No errors from deserialization failures
  expect_true(all(!res_shard$error))
})

test_that("shard backend parity holds under future::multisession", {
  skip_on_cran()
  skip_if_not_installed("future")

  ds <- gen_sample_dataset(c(5, 5, 5), 24, blocks = 3, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  # Baseline: sequential futures with shard backend
  baseline <- run_shard_with_plan(mspec, vox_iter, future::sequential)
  expect_equal(nrow(baseline), length(vox_iter))
  expect_true(all(!baseline$error))

  # Always test multisession
  res_multi <- run_shard_with_plan(mspec, vox_iter, future::multisession, workers = 2)
  expect_equal(res_multi$id, baseline$id)
  expect_equal(res_multi$error, baseline$error)
  if (!all(baseline$error)) {
    perf_base <- unlist(lapply(baseline$performance, function(p) if (is.null(p)) NA_real_ else p[[1]]))
    perf_multi <- unlist(lapply(res_multi$performance, function(p) if (is.null(p)) NA_real_ else p[[1]]))
    expect_equal(perf_multi, perf_base)
  }
})

test_that("shard backend parity holds under future::multicore", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not(isTRUE(future::supportsMulticore()), "future::multicore not supported in this environment")

  ds <- gen_sample_dataset(c(5, 5, 5), 24, blocks = 3, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  baseline <- run_shard_with_plan(mspec, vox_iter, future::sequential)
  res_mc <- run_shard_with_plan(mspec, vox_iter, future::multicore, workers = 2)

  expect_equal(res_mc$id, baseline$id)
  expect_equal(res_mc$error, baseline$error)
  if (!all(baseline$error)) {
    perf_base <- unlist(lapply(baseline$performance, function(p) if (is.null(p)) NA_real_ else p[[1]]))
    perf_mc <- unlist(lapply(res_mc$performance, function(p) if (is.null(p)) NA_real_ else p[[1]]))
    expect_equal(perf_mc, perf_base)
  }
})

test_that("shard backend parity holds under future.callr::callr", {
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.callr")

  ds <- gen_sample_dataset(c(5, 5, 5), 24, blocks = 3, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(5, length(vox_iter))]

  baseline <- run_shard_with_plan(mspec, vox_iter, future::sequential)
  res_callr <- run_shard_with_plan(mspec, vox_iter, future.callr::callr)

  expect_equal(res_callr$id, baseline$id)
  expect_equal(res_callr$error, baseline$error)
  if (!all(baseline$error)) {
    perf_base <- unlist(lapply(baseline$performance, function(p) if (is.null(p)) NA_real_ else p[[1]]))
    perf_callr <- unlist(lapply(res_callr$performance, function(p) if (is.null(p)) NA_real_ else p[[1]]))
    expect_equal(perf_callr, perf_base)
  }
})

test_that("switching future plans does not destabilize shard runs", {
  skip_on_cran()
  skip_if_not_installed("future")

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(4, length(vox_iter))]

  r1 <- run_shard_with_plan(mspec, vox_iter, future::sequential)
  r2 <- run_shard_with_plan(mspec, vox_iter, future::multisession, workers = 2)
  r3 <- run_shard_with_plan(mspec, vox_iter, future::sequential)

  expect_equal(r1$id, r2$id)
  expect_equal(r1$id, r3$id)
  expect_equal(r1$error, r2$error)
  expect_equal(r1$error, r3$error)
})

test_that("shard specs are single-use; reruns require fresh use_shard()", {
  skip_on_cran()
  skip_if_not_installed("future")

  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)

  sl <- get_searchlight(ds$dataset, radius = 3)
  vox_iter <- lapply(sl, function(x) x)
  vox_iter <- vox_iter[1:min(4, length(vox_iter))]

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)

  mspec_shard <- use_shard(mspec)

  set.seed(42)
  first_run <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )
  expect_true(all(!first_run$error))

  # Contract: mvpa_iterate() auto-cleans shard handles on exit, so a second
  # run with the same shard spec is not supported.
  set.seed(42)
  reuse_run <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard, vox_iter, ids = seq_along(vox_iter))
  )
  expect_true(all(reuse_run$error))
  expect_true(all(grepl("less than 2 features", reuse_run$error_message, fixed = TRUE)))

  # Supported path: refresh handles via use_shard() before rerunning.
  mspec_shard_fresh <- use_shard(mspec)
  set.seed(42)
  refreshed_run <- muffle_worker_version_warnings(
    mvpa_iterate(mspec_shard_fresh, vox_iter, ids = seq_along(vox_iter))
  )
  expect_true(all(!refreshed_run$error))
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

test_that("shard_extract_roi matches extract_roi for volumetric datasets", {
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
  mspec_shard <- use_shard(mspec)

  vox <- which(as.logical(mspec$dataset$mask))[1:6]
  sample <- get_samples(mspec$dataset, list(vox))$sample[[1]]

  roi_default <- rMVPA:::extract_roi(sample, mspec$dataset, min_voxels = 1)
  roi_shard <- rMVPA:::shard_extract_roi(vox, mspec_shard$shard_data, min_voxels = 1)

  expect_false(is.null(roi_default))
  expect_false(is.null(roi_shard))

  expect_equal(
    as.integer(neuroim2::indices(roi_default$train_roi)),
    as.integer(neuroim2::indices(roi_shard$train_roi))
  )
  expect_equal(
    neuroim2::values(roi_default$train_roi),
    neuroim2::values(roi_shard$train_roi)
  )

  expect_false(is.null(roi_default$test_roi))
  expect_false(is.null(roi_shard$test_roi))
  expect_equal(
    neuroim2::values(roi_default$test_roi),
    neuroim2::values(roi_shard$test_roi)
  )

  shard_cleanup(mspec_shard$shard_data)
})

test_that("shard_extract_roi matches extract_roi for clustered datasets", {
  ds <- gen_clustered_sample_dataset(
    D = c(8, 8, 8),
    nobs = 20,
    K = 8,
    nlevels = 2,
    blocks = 2,
    external_test = TRUE
  )
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec_shard <- use_shard(mspec)

  vox <- c(1L, 2L, 3L, 4L)
  sample <- get_samples(mspec$dataset, list(vox))$sample[[1]]

  roi_default <- rMVPA:::extract_roi(sample, mspec$dataset, min_voxels = 1)
  roi_shard <- rMVPA:::shard_extract_roi(vox, mspec_shard$shard_data, min_voxels = 1)

  expect_false(is.null(roi_default))
  expect_false(is.null(roi_shard))
  expect_equal(
    neuroim2::values(roi_default$train_roi),
    neuroim2::values(roi_shard$train_roi)
  )
  expect_equal(
    neuroim2::coords(roi_default$train_roi),
    neuroim2::coords(roi_shard$train_roi)
  )

  expect_false(is.null(roi_default$test_roi))
  expect_false(is.null(roi_shard$test_roi))
  expect_equal(
    neuroim2::values(roi_default$test_roi),
    neuroim2::values(roi_shard$test_roi)
  )

  shard_cleanup(mspec_shard$shard_data)
})

test_that("shard_extract_roi is equivalent to extract_roi across random volumetric neighborhoods", {
  set.seed(123)

  ds <- gen_sample_dataset(c(6, 6, 6), 24, blocks = 3, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec_shard <- use_shard(mspec)

  mask_idx <- which(as.logical(mspec$dataset$mask))
  max_idx <- prod(dim(mspec$dataset$mask)[1:3])
  invalid_pool <- c(0L, NA_integer_, max_idx + 1L, max_idx + 100L)

  for (i in seq_len(20)) {
    pool <- c(mask_idx, invalid_pool)
    nvox <- sample(1:20, 1)
    vox <- sample(pool, nvox, replace = TRUE)
    sample_obj <- get_samples(mspec$dataset, list(vox))$sample[[1]]

    roi_default <- rMVPA:::extract_roi(sample_obj, mspec$dataset, min_voxels = 1)
    roi_shard <- rMVPA:::shard_extract_roi(vox, mspec_shard$shard_data, min_voxels = 1)

    if (is.null(roi_default)) {
      expect_null(roi_shard)
    } else {
      expect_false(is.null(roi_shard))
      expect_equal(
        as.integer(neuroim2::indices(roi_default$train_roi)),
        as.integer(neuroim2::indices(roi_shard$train_roi))
      )
      expect_equal(
        neuroim2::values(roi_default$train_roi),
        neuroim2::values(roi_shard$train_roi)
      )
    }
  }

  shard_cleanup(mspec_shard$shard_data)
})

test_that("shard backend returns structured ROI errors for invalid neighborhoods", {
  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec <- use_shard(mspec)

  bad_vox <- c(0L, NA_integer_, prod(dim(ds$dataset$mask)) + 1000L)
  res <- muffle_worker_version_warnings(
    mvpa_iterate(mspec, list(bad_vox), ids = 1L)
  )

  expect_equal(nrow(res), 1L)
  expect_true(isTRUE(res$error[[1]]))
  expect_match(res$error_message[[1]], "failed validation|Error processing ROI")
  expect_true(isTRUE(res$warning[[1]]))
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

test_that("configure_runtime_backend removes shard state for backend='default'", {
  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec <- use_shard(mspec)

  mspec_default <- rMVPA:::configure_runtime_backend(mspec, backend = "default", context = "test")
  expect_false(inherits(mspec_default, "shard_model_spec"))
  expect_null(mspec_default$shard_data)

  mspec_default2 <- rMVPA:::configure_runtime_backend(mspec_default, backend = "default", context = "test")
  expect_false(inherits(mspec_default2, "shard_model_spec"))
  expect_null(mspec_default2$shard_data)
})

test_that("configure_runtime_backend backend='auto' falls back on unsupported datasets", {
  ds <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  mb <- mvpa_multibasis_dataset(
    train_data = list(ds$dataset$train_data, ds$dataset$train_data),
    mask = ds$dataset$mask
  )
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, mb, ds$design,
                      "classification", crossval = cval)

  resolved <- rMVPA:::configure_runtime_backend(mspec, backend = "auto", context = "test")
  expect_false(inherits(resolved, "shard_model_spec"))
  expect_null(resolved$shard_data)
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
    muffle_worker_version_warnings(
      mvpa_iterate(mspec, vox_iter, ids = seq_along(vox_iter))
    )
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
    muffle_worker_version_warnings(
      mvpa_iterate(mspec, vox_iter, ids = seq_along(vox_iter))
    )
  )
  expect_true(isTRUE(getOption("rMVPA.test.used_carrier_crate")))
})
