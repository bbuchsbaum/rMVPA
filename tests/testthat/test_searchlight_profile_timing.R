testthat::skip_if_not_installed("neuroim2")

build_profile_test_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 2, blocks = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
}

test_that("standard searchlight profiling attaches timing metadata", {
  mspec <- build_profile_test_spec()

  old_opt <- options(rMVPA.profile_searchlight = TRUE)
  on.exit(options(old_opt), add = TRUE)

  res <- run_searchlight(mspec, radius = 2, method = "standard")

  timing <- attr(res, "timing")
  expect_type(timing, "list")
  expect_true(is.numeric(timing$setup_seconds))
  expect_true(is.numeric(timing$iterate_seconds))
  expect_true(is.numeric(timing$combine_seconds))
  expect_true(is.list(timing$iterate))
  expect_true(is.list(timing$iterate$totals))
  expect_true(is.numeric(timing$iterate$totals$get_samples_seconds))
  expect_true(is.numeric(timing$iterate$totals$run_future_seconds))
})

test_that("profiling mode does not change searchlight output maps", {
  mspec <- build_profile_test_spec()

  set.seed(1201)
  old_off <- options(rMVPA.profile_searchlight = FALSE)
  on.exit(options(old_off), add = TRUE)
  res_base <- run_searchlight(mspec, radius = 2, method = "standard")

  set.seed(1201)
  old_on <- options(rMVPA.profile_searchlight = TRUE)
  on.exit(options(old_on), add = TRUE)
  res_profile <- run_searchlight(mspec, radius = 2, method = "standard")

  expect_searchlight_parity(
    reference = res_base,
    candidate = res_profile,
    atol = 1e-10,
    rtol = 1e-8
  )
})
