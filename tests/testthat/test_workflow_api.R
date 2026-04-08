library(rMVPA)

test_that("mvpa_config constructs public config objects", {
  cfg <- mvpa_config(
    mode = "searchlight",
    radius = 3,
    method = "standard"
  )

  expect_s3_class(cfg, "rmvpa_config")
  expect_identical(cfg$mode, "searchlight")
  expect_identical(cfg$radius, 3)
})

test_that("build_analysis accepts a prebuilt model spec", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 2)
  mspec <- mvpa_model(
    load_model("corclass"),
    dataset = ds$dataset,
    design = ds$design,
    crossval = blocked_cross_validation(ds$design$block_var)
  )

  analysis <- build_analysis(mvpa_config(
    mode = "searchlight",
    model_spec = mspec,
    radius = 2
  ))

  expect_s3_class(analysis, "rmvpa_analysis")
  expect_length(analysis$entries, 1)
  expect_identical(analysis$entries[[1]]$model_spec, mspec)
})

test_that("run_analysis executes configured workflows and attaches context", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 3, blocks = 2)
  mspec <- mvpa_model(
    load_model("corclass"),
    dataset = ds$dataset,
    design = ds$design,
    crossval = blocked_cross_validation(ds$design$block_var)
  )

  analysis <- build_analysis(mvpa_config(
    mode = "searchlight",
    model_spec = mspec,
    radius = 2
  ))

  out <- testthat::with_mocked_bindings(
    run_analysis(analysis, preflight = "warn"),
    run_searchlight = function(model_spec, radius, method, niter, preflight, ...) {
      structure(
        list(
          results = list(),
          n_voxels = 0L,
          active_voxels = 0L,
          metrics = character(0)
        ),
        class = c("searchlight_result", "list")
      )
    },
    .package = "rMVPA"
  )

  expect_s3_class(out, "rmvpa_analysis_run")
  expect_length(out$results, 1)
  expect_equal(attr(out$results[[1]], "analysis_context")$mode, "searchlight")
})

test_that("save_results handles rmvpa_analysis_run objects", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 3, blocks = 2)
  mspec <- mvpa_model(
    load_model("corclass"),
    dataset = ds$dataset,
    design = ds$design,
    crossval = blocked_cross_validation(ds$design$block_var)
  )

  analysis <- build_analysis(mvpa_config(
    mode = "searchlight",
    model_spec = mspec,
    radius = 2
  ))

  run_obj <- structure(
    list(
      analysis = analysis,
      results = list(
        analysis = structure(
          list(
            results = list(),
            n_voxels = 0L,
            active_voxels = 0L,
            metrics = character(0)
          ),
          class = c("searchlight_result", "list")
        )
      ),
      preflight = "warn"
    ),
    class = c("rmvpa_analysis_run", "list")
  )

  tmp <- tempfile()
  save_results(run_obj, dir = tmp, quiet = TRUE)

  expect_true(file.exists(file.path(tmp, "aux", "analysis_config.rds")))
  expect_true(length(list.files(tmp, pattern = "^manifest\\.")) > 0L)

  unlink(tmp, recursive = TRUE)
})
