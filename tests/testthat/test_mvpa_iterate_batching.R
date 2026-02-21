testthat::skip_if_not_installed("neuroim2")

build_iterate_test_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 18, nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
  list(mspec = mspec, dataset = ds)
}

make_passthrough_processor <- function() {
  function(obj, roi, rnum, center_global_id = NA) {
    tibble::tibble(
      result = list(NULL),
      indices = list(NULL),
      performance = list(NULL),
      id = rnum,
      error = FALSE,
      error_message = "~",
      warning = FALSE,
      warning_message = "~"
    )
  }
}

test_that("mvpa_iterate batch profiling keeps processed/skipped accounting stable", {
  built <- build_iterate_test_spec()
  mspec <- built$mspec
  mask_idx <- which(built$dataset$dataset$mask > 0)
  vox_list <- list(
    as.integer(mask_idx[1:3]),
    as.integer(mask_idx[4:6]),
    integer(0),
    integer(0)
  )
  ids <- c(101L, 102L, 103L, 104L)

  old_opt <- options(
    rMVPA.profile_searchlight = TRUE,
    rMVPA.searchlight_backend_default = "default"
  )
  on.exit(options(old_opt), add = TRUE)

  res <- mvpa_iterate(
    mod_spec = mspec,
    vox_list = vox_list,
    ids = ids,
    batch_size = 2,
    verbose = FALSE,
    analysis_type = "searchlight",
    processor = make_passthrough_processor(),
    fail_fast = FALSE
  )

  timing <- attr(res, "timing")
  expect_type(timing, "list")
  expect_identical(as.integer(timing$processed_rois), 2L)
  expect_identical(as.integer(timing$skipped_rois), 2L)
  expect_identical(nrow(res), 4L)

  skipped_rows <- res[res$id %in% c(103L, 104L), , drop = FALSE]
  expect_identical(nrow(skipped_rows), 2L)
  expect_true(all(skipped_rows$error))
  expect_true(all(grepl("ROI filtered out", skipped_rows$error_message, fixed = TRUE)))
})

test_that("mvpa_iterate is stable under fail_fast toggle when processor succeeds", {
  built <- build_iterate_test_spec()
  mspec <- built$mspec
  mask_idx <- as.integer(which(built$dataset$dataset$mask > 0))

  vox_list <- list(
    mask_idx[1:4],
    mask_idx[5:8],
    mask_idx[9:12],
    mask_idx[13:16]
  )
  ids <- c(201L, 202L, 203L, 204L)
  processor <- make_passthrough_processor()

  old_opt <- options(rMVPA.searchlight_backend_default = "default")
  on.exit(options(old_opt), add = TRUE)

  res_slow <- mvpa_iterate(
    mod_spec = mspec,
    vox_list = vox_list,
    ids = ids,
    batch_size = 2,
    verbose = FALSE,
    analysis_type = "searchlight",
    processor = processor,
    fail_fast = FALSE
  )

  res_fast <- mvpa_iterate(
    mod_spec = mspec,
    vox_list = vox_list,
    ids = ids,
    batch_size = 2,
    verbose = FALSE,
    analysis_type = "searchlight",
    processor = processor,
    fail_fast = TRUE
  )

  expect_equal(res_slow$id, res_fast$id)
  expect_equal(res_slow$error, res_fast$error)
  expect_equal(res_slow$warning, res_fast$warning)
  expect_equal(res_slow$error_message, res_fast$error_message)
  expect_equal(res_slow$warning_message, res_fast$warning_message)
})
