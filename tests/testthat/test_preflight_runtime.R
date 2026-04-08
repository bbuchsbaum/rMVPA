library(rMVPA)

test_that("run_searchlight.default stops when preflight policy is error", {
  expect_error(
    testthat::with_mocked_bindings(
      rMVPA:::run_searchlight.default(
        model_spec = structure(list(), class = "mvpa_model"),
        radius = 2,
        method = "standard",
        preflight = "error"
      ),
      .plugin_preflight = function(...) NULL,
      .apply_analysis_preflight = function(...) {
        stop("synthetic preflight failure", call. = FALSE)
      },
      .run_searchlight_engine = function(...) {
        structure(
          list(handled = TRUE, result = list(), engine = "swift"),
          class = "list"
        )
      },
      .package = "rMVPA"
    ),
    "synthetic preflight failure"
  )
})

test_that("apply_analysis_preflight returns validation results for supported specs", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 5)
  mspec <- mvpa_model(
    load_model("corclass"),
    dataset = ds$dataset,
    design = ds$design,
    crossval = blocked_cross_validation(ds$design$block_var)
  )

  out <- rMVPA:::.apply_analysis_preflight(mspec, preflight = "warn", context = "unit-test")
  expect_s3_class(out, "validation_result")
})
