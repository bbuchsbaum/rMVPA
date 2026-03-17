test_that("create_model_spec is exported and validates core inputs", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 12)

  mspec <- create_model_spec("toy_model_spec", ds$dataset, ds$design)
  expect_s3_class(mspec, "model_spec")
  expect_true("toy_model_spec" %in% class(mspec))

  expect_error(
    create_model_spec("", ds$dataset, ds$design),
    "name"
  )
  expect_error(
    create_model_spec("toy_model_spec", list(), ds$design),
    "mvpa_dataset"
  )
  expect_error(
    create_model_spec("toy_model_spec", ds$dataset, list()),
    "design"
  )
})

test_that("y_train/y_test model_spec methods fall back to design", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 16, external_test = TRUE)
  mspec <- create_model_spec("toy_model_spec", ds$dataset, ds$design)

  expect_equal(y_train(mspec), y_train(ds$design))
  expect_equal(y_test(mspec), y_test(ds$design))
})

test_that("minimal fit_roi plugin works without custom y_train method", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 20, blocks = 2, nlevels = 2)
  cv <- blocked_cross_validation(ds$design$block_var)

  mspec <- create_model_spec(
    "toy_plugin_model",
    dataset = ds$dataset,
    design = ds$design,
    crossval = cv,
    compute_performance = TRUE,
    return_predictions = FALSE
  )

  fit_roi.toy_plugin_model <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(score = mean(roi_data$train_data)),
      indices = roi_data$indices,
      id = context$id
    )
  }
  output_schema.toy_plugin_model <- function(model) {
    list(score = "scalar")
  }

  registerS3method(
    "fit_roi",
    "toy_plugin_model",
    fit_roi.toy_plugin_model,
    envir = asNamespace("rMVPA")
  )
  registerS3method(
    "output_schema",
    "toy_plugin_model",
    output_schema.toy_plugin_model,
    envir = asNamespace("rMVPA")
  )

  out <- suppressWarnings(run_searchlight(mspec, radius = 2, method = "standard"))
  expect_s3_class(out, "searchlight_result")
  expect_true("score" %in% out$metrics)
})
