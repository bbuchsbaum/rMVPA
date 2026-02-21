testthat::skip_if_not_installed("neuroim2")

test_that("geometry cache is enabled by default", {
  expect_true(rMVPA:::.searchlight_geometry_cache_enabled())
})

build_geometry_cache_dataset <- function() {
  ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 16, nlevels = 2, blocks = 4)
  ds$dataset
}

window_signature <- function(win) {
  coords <- as.matrix(win@coords)
  idx <- seq_len(nrow(coords))
  sig <- sum(
    coords[, 1] * (idx + 1) +
      coords[, 2] * (idx + 7) * 3 +
      coords[, 3] * (idx + 13) * 5
  )
  paste(
    as.integer(win@parent_index),
    as.integer(win@center_index),
    nrow(coords),
    formatC(sig, format = "f", digits = 0),
    sep = ":"
  )
}

test_that("standard image searchlight geometry cache reuses computed neighborhoods", {
  dset <- build_geometry_cache_dataset()
  base_compute <- rMVPA:::.compute_image_searchlight
  call_count <- 0L

  rMVPA:::.searchlight_geometry_cache_clear()

  testthat::with_mocked_bindings(
    {
      sl1 <- get_searchlight(dset, type = "standard", radius = 2)
      sl2 <- get_searchlight(dset, type = "standard", radius = 2)
      expect_equal(length(sl1), length(sl2))
    },
    .compute_image_searchlight = function(...) {
      call_count <<- call_count + 1L
      base_compute(...)
    },
    .package = "rMVPA"
  )

  expect_identical(call_count, 1L)
  expect_identical(rMVPA:::.searchlight_geometry_cache_size(), 1L)
})

test_that("standard image searchlight geometry cache is always active", {
  dset <- build_geometry_cache_dataset()
  base_compute <- rMVPA:::.compute_image_searchlight
  call_count <- 0L
  rMVPA:::.searchlight_geometry_cache_clear()

  testthat::with_mocked_bindings(
    {
      get_searchlight(dset, type = "standard", radius = 2)
      get_searchlight(dset, type = "standard", radius = 2)
    },
    .compute_image_searchlight = function(...) {
      call_count <<- call_count + 1L
      base_compute(...)
    },
    .package = "rMVPA"
  )

  expect_identical(call_count, 1L)
})

test_that("geometry cache invalidates on mask changes", {
  dset1 <- build_geometry_cache_dataset()
  mask_vals <- as.numeric(neuroim2::values(dset1$mask))
  mask_vals[1] <- 0
  mask_arr <- array(mask_vals > 0, dim = dim(dset1$mask)[1:3])
  mask2 <- neuroim2::LogicalNeuroVol(mask_arr, neuroim2::space(dset1$mask))
  dset2 <- mvpa_dataset(train_data = dset1$train_data, mask = mask2)

  base_compute <- rMVPA:::.compute_image_searchlight
  call_count <- 0L

  rMVPA:::.searchlight_geometry_cache_clear()

  testthat::with_mocked_bindings(
    {
      get_searchlight(dset1, type = "standard", radius = 2)
      get_searchlight(dset2, type = "standard", radius = 2)
    },
    .compute_image_searchlight = function(...) {
      call_count <<- call_count + 1L
      base_compute(...)
    },
    .package = "rMVPA"
  )

  expect_identical(call_count, 2L)
})

test_that("cached geometry path preserves standard searchlight membership", {
  dset <- build_geometry_cache_dataset()

  rMVPA:::.searchlight_geometry_cache_clear()
  sl_base <- get_searchlight(dset, type = "standard", radius = 2)
  rMVPA:::.searchlight_geometry_cache_clear()
  sl_cached <- get_searchlight(dset, type = "standard", radius = 2)

  expect_equal(length(sl_base), length(sl_cached))
  sig_base <- vapply(sl_base, window_signature, FUN.VALUE = character(1))
  sig_cached <- vapply(sl_cached, window_signature, FUN.VALUE = character(1))
  expect_equal(sig_base, sig_cached)
})

test_that("non-standard searchlights bypass geometry cache", {
  dset <- build_geometry_cache_dataset()
  base_compute <- rMVPA:::.compute_image_searchlight
  call_count <- 0L
  rMVPA:::.searchlight_geometry_cache_clear()

  testthat::with_mocked_bindings(
    {
      get_searchlight(dset, type = "resampled", radius = 2, iter = 1L)
      get_searchlight(dset, type = "resampled", radius = 2, iter = 1L)
    },
    .compute_image_searchlight = function(...) {
      call_count <<- call_count + 1L
      base_compute(...)
    },
    .package = "rMVPA"
  )

  expect_identical(call_count, 2L)
  expect_identical(rMVPA:::.searchlight_geometry_cache_size(), 0L)
})

test_that("geometry cache retains multiple keys under default LRU budget", {
  dset <- build_geometry_cache_dataset()
  base_compute <- rMVPA:::.compute_image_searchlight
  call_count <- 0L
  rMVPA:::.searchlight_geometry_cache_clear()

  testthat::with_mocked_bindings(
    {
      get_searchlight(dset, type = "standard", radius = 2)
      get_searchlight(dset, type = "standard", radius = 3)
      get_searchlight(dset, type = "standard", radius = 2)
    },
    .compute_image_searchlight = function(...) {
      call_count <<- call_count + 1L
      base_compute(...)
    },
    .package = "rMVPA"
  )

  expect_identical(call_count, 2L)
  expect_identical(rMVPA:::.searchlight_geometry_cache_size(), 2L)
})
