test_that("mock_roi_data builds expected matrix payload", {
  roi_data <- mock_roi_data(n_train = 12, n_features = 5, n_test = 4, seed = 123)

  expect_true(is.matrix(roi_data$train_data))
  expect_equal(dim(roi_data$train_data), c(12, 5))
  expect_true(is.matrix(roi_data$test_data))
  expect_equal(dim(roi_data$test_data), c(4, 5))
  expect_equal(roi_data$indices, as.integer(1:5))
  expect_null(roi_data$train_roi)
  expect_null(roi_data$test_roi)
})

test_that("mock_context provides a valid default fit_roi context", {
  ctx <- mock_context(id = 11L)

  expect_true(inherits(ctx$design, "mvpa_design"))
  expect_true(is.null(ctx$cv_spec))
  expect_equal(ctx$id, 11L)
  expect_true(is.na(ctx$center_global_id))
})

test_that("validate_plugin_model accepts a valid toy plugin contract", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 16, nlevels = 2)
  model_spec <- create_model_spec("toy_validate_plugin_model", ds$dataset, ds$design)

  fit_roi.toy_validate_plugin_model <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(score = mean(roi_data$train_data)),
      indices = roi_data$indices,
      id = context$id
    )
  }

  output_schema.toy_validate_plugin_model <- function(model) {
    list(score = "scalar")
  }

  registerS3method(
    "fit_roi",
    "toy_validate_plugin_model",
    fit_roi.toy_validate_plugin_model,
    envir = asNamespace("rMVPA")
  )
  registerS3method(
    "output_schema",
    "toy_validate_plugin_model",
    output_schema.toy_validate_plugin_model,
    envir = asNamespace("rMVPA")
  )

  out <- validate_plugin_model(
    model_spec,
    roi_data = mock_roi_data(n_train = 16, n_features = 6, seed = 1),
    context = mock_context(design = ds$design, id = 7L)
  )

  expect_s3_class(out, "plugin_validation_result")
  expect_true(out$valid)
  expect_equal(out$model_class, "toy_validate_plugin_model")
  expect_equal(out$metric_names, "score")
  expect_equal(out$metric_count, 1L)
})

test_that("validate_plugin_model errors on schema/metric mismatch", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 16, nlevels = 2)
  model_spec <- create_model_spec("toy_validate_plugin_mismatch", ds$dataset, ds$design)

  fit_roi.toy_validate_plugin_mismatch <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(other_score = mean(roi_data$train_data)),
      indices = roi_data$indices,
      id = context$id
    )
  }

  output_schema.toy_validate_plugin_mismatch <- function(model) {
    list(score = "scalar")
  }

  registerS3method(
    "fit_roi",
    "toy_validate_plugin_mismatch",
    fit_roi.toy_validate_plugin_mismatch,
    envir = asNamespace("rMVPA")
  )
  registerS3method(
    "output_schema",
    "toy_validate_plugin_mismatch",
    output_schema.toy_validate_plugin_mismatch,
    envir = asNamespace("rMVPA")
  )

  expect_error(
    validate_plugin_model(
      model_spec,
      roi_data = mock_roi_data(n_train = 16, n_features = 6, seed = 1),
      context = mock_context(design = ds$design, id = 3L)
    ),
    "name mismatch"
  )
})

test_that("cv_evaluate_roi runs internal CV and can return roi_result or row", {
  ds <- gen_sample_dataset(c(5, 5, 5), nobs = 24, blocks = 3, nlevels = 2)
  cv <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification", crossval = cv)

  vox <- sample(which(ds$dataset$mask > 0), 25)
  roi <- as_roi(data_sample(ds$dataset, vox), ds$dataset)
  roi_data <- list(
    train_data = as.matrix(neuroim2::values(roi$train_roi)),
    test_data = NULL,
    indices = neuroim2::indices(roi$train_roi),
    train_roi = roi$train_roi,
    test_roi = roi$test_roi
  )
  context <- mock_context(design = ds$design, cv_spec = cv, id = 9L)

  out <- cv_evaluate_roi(mspec, roi_data, context, mode = "internal")
  expect_s3_class(out, "roi_result")
  expect_equal(out$id, 9L)

  row <- cv_evaluate_roi(mspec, roi_data, context, mode = "internal", return = "row")
  expect_s3_class(row, "tbl_df")
  expect_true(all(c("result", "performance", "id", "error") %in% names(row)))
})

test_that("cv_evaluate_roi external mode requires model spec with external labels", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 20, nlevels = 2)
  mspec <- create_model_spec("toy_cv_external_guard", ds$dataset, ds$design)
  roi_data <- mock_roi_data(n_train = 20, n_test = 5, n_features = 4, seed = 42)

  expect_error(
    cv_evaluate_roi(
      mspec,
      roi_data = roi_data,
      context = mock_context(design = ds$design, id = 1L),
      mode = "external"
    ),
    "requires a model spec with external test labels"
  )
})

test_that("validate_model_spec reports pass/fail states", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 16, nlevels = 2)

  good <- create_model_spec("toy_validate_model_spec", ds$dataset, ds$design)
  fit_roi.toy_validate_model_spec <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(score = mean(roi_data$train_data)),
      indices = roi_data$indices,
      id = context$id
    )
  }
  output_schema.toy_validate_model_spec <- function(model) {
    list(score = "scalar")
  }
  registerS3method("fit_roi", "toy_validate_model_spec",
                   fit_roi.toy_validate_model_spec, envir = asNamespace("rMVPA"))
  registerS3method("output_schema", "toy_validate_model_spec",
                   output_schema.toy_validate_model_spec, envir = asNamespace("rMVPA"))

  vgood <- validate_model_spec(
    good,
    require_schema = TRUE,
    dry_run = TRUE,
    roi_data = mock_roi_data(n_train = 16, n_features = 4, seed = 2),
    context = mock_context(design = ds$design, id = 1L)
  )
  expect_s3_class(vgood, "model_spec_validation_result")
  expect_true(vgood$valid)

  bad <- create_model_spec("toy_validate_model_spec_missing_fit", ds$dataset, ds$design)
  vbad <- validate_model_spec(bad, require_schema = FALSE, dry_run = FALSE)
  expect_s3_class(vbad, "model_spec_validation_result")
  expect_false(vbad$valid)
  statuses <- vapply(vbad$checks, `[[`, character(1), "status")
  expect_true(any(statuses == "fail"))
})

test_that("validate_model_spec resolves fit_roi methods through the class chain", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 16, nlevels = 2)
  mspec <- create_model_spec("toy_parent_validate_model_spec", ds$dataset, ds$design)
  class(mspec) <- c("toy_child_validate_model_spec", class(mspec))

  fit_roi.toy_parent_validate_model_spec <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(score = mean(roi_data$train_data)),
      indices = roi_data$indices,
      id = context$id
    )
  }
  registerS3method(
    "fit_roi",
    "toy_parent_validate_model_spec",
    fit_roi.toy_parent_validate_model_spec,
    envir = asNamespace("rMVPA")
  )

  out <- validate_model_spec(mspec, require_schema = FALSE, dry_run = FALSE)
  expect_s3_class(out, "model_spec_validation_result")
  expect_true(out$valid)

  check_names <- vapply(out$checks, `[[`, character(1), "name")
  check_msgs <- vapply(out$checks, `[[`, character(1), "message")
  fit_idx <- which(check_names == "fit_roi_method")
  expect_true(length(fit_idx) == 1L)
  expect_match(check_msgs[fit_idx], "fit_roi\\.toy_parent_validate_model_spec")
})

test_that("run_searchlight fails fast on plugin preflight failure", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 16, nlevels = 2)
  mspec <- create_model_spec("toy_preflight_missing_fit_search", ds$dataset, ds$design)

  expect_error(
    run_searchlight(mspec, radius = 2, method = "standard"),
    "plugin contract preflight failed"
  )
})

test_that("run_regional fails fast on plugin preflight failure", {
  ds <- gen_sample_dataset(c(4, 4, 4), nobs = 16, nlevels = 2)
  mspec <- create_model_spec("toy_preflight_missing_fit_regional", ds$dataset, ds$design)

  region_mask <- NeuroVol(
    sample(1:2, size = length(ds$dataset$mask), replace = TRUE),
    space(ds$dataset$mask)
  )

  expect_error(
    run_regional(mspec, region_mask),
    "plugin contract preflight failed"
  )
})
