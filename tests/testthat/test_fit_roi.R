test_that("roi_result constructs correctly", {
  res <- roi_result(
    metrics = c(accuracy = 0.85, AUC = 0.9),
    indices = 1:10,
    id = 42
  )
  expect_s3_class(res, "roi_result")
  expect_equal(res$metrics, c(accuracy = 0.85, AUC = 0.9))
  expect_equal(res$indices, 1:10)
  expect_equal(res$id, 42)
  expect_null(res$result)
  expect_false(res$error)
  expect_equal(res$error_message, "~")
})

test_that("roi_result handles error case", {
  res <- roi_result(
    metrics = NULL,
    indices = 1:5,
    id = 7,
    error = TRUE,
    error_message = "Too few features"
  )
  expect_s3_class(res, "roi_result")
  expect_true(res$error)
  expect_equal(res$error_message, "Too few features")
  expect_null(res$metrics)
})

test_that("roi_result_to_tibble converts success", {
  res <- roi_result(
    metrics = c(accuracy = 0.85),
    indices = 1:10,
    id = 42,
    result = list(some_detail = "ok")
  )
  tbl <- rMVPA:::roi_result_to_tibble(res)
  expect_s3_class(tbl, "tbl_df")
  expect_equal(nrow(tbl), 1)
  expect_equal(tbl$id, 42)
  expect_false(tbl$error)
  expect_equal(tbl$error_message, "~")
  expect_equal(tbl$performance[[1]], c(accuracy = 0.85))
  expect_equal(tbl$indices[[1]], 1:10)
  expect_equal(tbl$result[[1]], list(some_detail = "ok"))
})

test_that("roi_result_to_tibble converts error", {
  res <- roi_result(
    metrics = NULL,
    indices = 1:5,
    id = 7,
    error = TRUE,
    error_message = "boom"
  )
  tbl <- rMVPA:::roi_result_to_tibble(res)
  expect_equal(nrow(tbl), 1)
  expect_true(tbl$error)
  expect_equal(tbl$error_message, "boom")
  expect_null(tbl$performance[[1]])
})

test_that("output_schema.mvpa_model returns scalar schema for standard model", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
  schema <- output_schema(mspec)
  expect_true(is.list(schema))
  expect_equal(names(schema), c("Accuracy", "AUC"))
  expect_true(all(unlist(schema) == "scalar"))
})

test_that("fit_roi dispatches to model-specific method", {
  # Define a toy model class
  toy <- structure(list(), class = c("toy_fit_roi_model", "list"))

  # Register a fit_roi method for the toy class
  fit_roi.toy_fit_roi_model <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(toy_metric = 1.0),
      indices = roi_data$indices,
      id = context$id,
      result = list(note = "toy fit")
    )
  }

  registerS3method("fit_roi", "toy_fit_roi_model", fit_roi.toy_fit_roi_model,
                   envir = asNamespace("rMVPA"))

  # Verify .has_fit_roi detects the method
  expect_true(rMVPA:::.has_fit_roi(toy))

  # Call fit_roi and verify result
  roi_data <- list(
    train_data = matrix(rnorm(20), 4, 5),
    test_data = NULL,
    indices = 1:5,
    train_roi = NULL,
    test_roi = NULL
  )
  context <- list(design = NULL, cv_spec = NULL, id = 99, center_global_id = NA)

  res <- fit_roi(toy, roi_data, context)
  expect_s3_class(res, "roi_result")
  expect_equal(res$metrics, c(toy_metric = 1.0))
  expect_equal(res$id, 99)
  expect_false(res$error)
})

test_that(".has_fit_roi returns TRUE for mvpa_model (Phase 4)", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 2)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = blocked_cross_validation(ds$design$block_var))
  expect_true(rMVPA:::.has_fit_roi(mspec))
})

test_that("process_roi.default still works for existing models (regression test)", {
  # Use a larger dataset to avoid classifier edge cases with small ROIs
  ds <- gen_sample_dataset(c(6, 6, 6), 80, blocks = 4, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = cval)

  # Get voxels and construct an ROI
  vox <- sample(which(ds$dataset$mask > 0), 30)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)

  # process_roi should use the legacy path (internal_crossval)
  result <- process_roi(mspec, roi, 1)
  expect_s3_class(result, "tbl_df")
  expect_true("error" %in% names(result))
  expect_true("performance" %in% names(result))
  expect_true("id" %in% names(result))
})

test_that("process_roi.default prefers fit_roi when method exists", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = cval)

  # Add a custom class that has a fit_roi method
  class(mspec) <- c("test_shim_model", class(mspec))

  fit_roi.test_shim_model <- function(model, roi_data, context, ...) {
    roi_result(
      metrics = c(shim_test = 42.0),
      indices = roi_data$indices,
      id = context$id
    )
  }
  registerS3method("fit_roi", "test_shim_model", fit_roi.test_shim_model,
                   envir = asNamespace("rMVPA"))

  vox <- sample(which(ds$dataset$mask > 0), 15)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)

  result <- process_roi(mspec, roi, 1)
  expect_s3_class(result, "tbl_df")
  expect_false(result$error)
  # The performance should come from our fit_roi method
  expect_equal(result$performance[[1]], c(shim_test = 42.0))
})

test_that("process_roi.default handles fit_roi errors gracefully", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = cval)

  class(mspec) <- c("test_error_model", class(mspec))

  fit_roi.test_error_model <- function(model, roi_data, context, ...) {
    stop("intentional test error")
  }
  registerS3method("fit_roi", "test_error_model", fit_roi.test_error_model,
                   envir = asNamespace("rMVPA"))

  vox <- sample(which(ds$dataset$mask > 0), 15)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)

  result <- process_roi(mspec, roi, 1)
  expect_s3_class(result, "tbl_df")
  expect_true(result$error)
  expect_match(result$error_message, "intentional test error")
})

# ---- Phase 4 Tests ----

test_that("generate_folds delegates to crossval_samples", {
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  y <- ds$design$y_train
  data <- tibble::tibble(.row = seq_along(y))

  folds <- generate_folds(cval, data, y)
  expect_s3_class(folds, "tbl_df")
  expect_true(all(c("ytrain", "ytest", "train", "test", ".id") %in% names(folds)))
  expect_true(nrow(folds) > 0)
})

test_that("fit_roi.mvpa_model returns roi_result", {
  ds <- gen_sample_dataset(c(6, 6, 6), 80, blocks = 4, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = cval)

  vox <- sample(which(ds$dataset$mask > 0), 30)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)

  roi_data <- list(
    train_data = NULL,
    test_data = NULL,
    indices = vox,
    train_roi = roi$train_roi,
    test_roi = roi$test_roi
  )
  context <- list(design = mspec$design, cv_spec = mspec$crossval,
                  id = 1, center_global_id = NA)

  res <- fit_roi(mspec, roi_data, context)
  expect_s3_class(res, "roi_result")
  expect_false(res$error)
  expect_equal(res$id, 1)
  expect_equal(res$indices, vox)
  expect_true(!is.null(res$metrics))
})

test_that("fit_roi.mvpa_model via process_roi.default produces valid tibble", {
  ds <- gen_sample_dataset(c(6, 6, 6), 80, blocks = 4, nlevels = 2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = cval)

  vox <- sample(which(ds$dataset$mask > 0), 30)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)

  result <- process_roi(mspec, roi, 1)
  expect_s3_class(result, "tbl_df")
  expect_true("error" %in% names(result))
  expect_true("performance" %in% names(result))
  expect_true("id" %in% names(result))
  expect_false(result$error)
})
