test_that("comp_perf errors on regional schema/performance width mismatch", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  model_spec <- structure(
    list(dataset = ds$dataset),
    class = c("toy_regional_schema_width_model", "model_spec", "list")
  )

  output_schema.toy_regional_schema_width_model <- function(model) {
    list(metric_a = "scalar", metric_b = "scalar")
  }
  registerS3method(
    "output_schema",
    "toy_regional_schema_width_model",
    output_schema.toy_regional_schema_width_model,
    envir = asNamespace("rMVPA")
  )

  region_mask <- NeuroVol(
    sample(1:2, size = length(ds$dataset$mask), replace = TRUE),
    space(ds$dataset$mask)
  )
  results <- tibble::tibble(
    id = c(1L, 2L),
    performance = list(c(metric_a = 0.1), c(metric_a = 0.2))
  )

  expect_error(
    rMVPA:::comp_perf(results, region_mask, model_spec = model_spec),
    "width mismatch"
  )
})

test_that("comp_perf errors on regional schema/performance name mismatch", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  model_spec <- structure(
    list(dataset = ds$dataset),
    class = c("toy_regional_schema_name_model", "model_spec", "list")
  )

  output_schema.toy_regional_schema_name_model <- function(model) {
    list(metric_a = "scalar")
  }
  registerS3method(
    "output_schema",
    "toy_regional_schema_name_model",
    output_schema.toy_regional_schema_name_model,
    envir = asNamespace("rMVPA")
  )

  region_mask <- NeuroVol(
    sample(1:2, size = length(ds$dataset$mask), replace = TRUE),
    space(ds$dataset$mask)
  )
  results <- tibble::tibble(
    id = c(1L, 2L),
    performance = list(c(other_name = 0.1), c(other_name = 0.2))
  )

  expect_error(
    rMVPA:::comp_perf(results, region_mask, model_spec = model_spec),
    "metric name mismatch"
  )
})

test_that("run_regional degrades gracefully when all ROI schema checks fail", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 18, nlevels = 2)
  model_spec <- create_model_spec(
    "toy_regional_bad_schema_model",
    dataset = ds$dataset,
    design = ds$design,
    compute_performance = TRUE,
    return_predictions = FALSE
  )

  fit_roi.toy_regional_bad_schema_model <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(other_metric = mean(roi_data$train_data)),
      indices = roi_data$indices,
      id = context$id
    )
  }
  output_schema.toy_regional_bad_schema_model <- function(model) {
    list(score = "scalar")
  }

  registerS3method(
    "fit_roi",
    "toy_regional_bad_schema_model",
    fit_roi.toy_regional_bad_schema_model,
    envir = asNamespace("rMVPA")
  )
  registerS3method(
    "output_schema",
    "toy_regional_bad_schema_model",
    output_schema.toy_regional_bad_schema_model,
    envir = asNamespace("rMVPA")
  )

  region_mask <- NeuroVol(
    sample(1:2, size = length(ds$dataset$mask), replace = TRUE),
    space(ds$dataset$mask)
  )

  res <- run_regional(model_spec, region_mask)
  expect_s3_class(res, "regional_mvpa_result")
  expect_true(is.data.frame(res$performance_table))
  expect_equal(names(res$performance_table), "roinum")
  expect_equal(length(res$vol_results), 0L)
})
