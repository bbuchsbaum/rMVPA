context("remap_rrr_model")

test_that("remap_rrr_model runs with identity fallback (rank=0)", {
  skip_on_cran()
  set.seed(1)
  toy <- gen_sample_dataset(D = c(5,5,5), nobs = 60, nlevels = 3, blocks = 3, external_test = TRUE)
  # simple regional mask: 3 ROIs randomly
  regionMask <- neuroim2::NeuroVol(sample(1:3, size = length(toy$dataset$mask), replace = TRUE),
                                   neuroim2::space(toy$dataset$mask))

  mspec <- remap_rrr_model(
    dataset = toy$dataset,
    design  = toy$design,
    base_classifier = "sda_notune",
    link_by = NULL,
    rank = 0,  # force identity/no-adapt for fast, dependency-free run
    leave_one_key_out = TRUE
  )

  res <- run_regional(mspec, regionMask)
  expect_s3_class(res, "regional_mvpa_result")
  # basic fields present
  expect_true(is.list(res$vol_results))
})

test_that("remap_rrr_model LOKO path executes when rank>0 but rrpack missing", {
  skip_on_cran()
  set.seed(2)
  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 48, nlevels = 3, blocks = 3, external_test = TRUE)
  mask <- neuroim2::NeuroVol(sample(1:2, size = length(toy$dataset$mask), replace = TRUE),
                             neuroim2::space(toy$dataset$mask))
  # Ask for auto rank but if rrpack is unavailable, code should fallback gracefully
  mspec <- remap_rrr_model(
    dataset = toy$dataset,
    design  = toy$design,
    base_classifier = "sda_notune",
    rank = "auto",
    max_rank = 3,
    leave_one_key_out = TRUE
  )
  res <- run_regional(mspec, mask)
  expect_s3_class(res, "regional_mvpa_result")
})

test_that("remap_rrr_model runs randomized searchlight with niter=4", {
  skip_on_cran()
  set.seed(313)
  toy <- gen_sample_dataset(D = c(4,4,4), nobs = 48, nlevels = 3, blocks = 3, external_test = TRUE)

  mspec <- remap_rrr_model(
    dataset = toy$dataset,
    design  = toy$design,
    rank = 0,                 # avoid rrpack dependency; just test searchlight plumbing
    leave_one_key_out = FALSE # faster for randomized searchlight
  )

  res <- run_searchlight(
    mspec,
    radius = 2,
    method = "randomized",
    niter = 4
  )

  expect_s3_class(res, "searchlight_result")
  expect_gt(length(res$results), 0)
  expect_true(any(!is.na(neuroim2::values(res$results[[1]]))))
})
