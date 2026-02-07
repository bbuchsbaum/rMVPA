test_that("as_worker_spec removes dataset and preserves class chain", {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")

  ws <- rMVPA:::as_worker_spec(mspec)

  expect_null(ws$dataset)
  expect_equal(class(ws), class(mspec))
  expect_identical(ws$design, mspec$design)
  expect_identical(ws$model, mspec$model)
  expect_identical(ws$crossval, mspec$crossval)
  expect_identical(ws$performance, mspec$performance)
  expect_identical(ws$return_predictions, mspec$return_predictions)
  expect_identical(ws$compute_performance, mspec$compute_performance)
})

test_that("as_worker_spec is idempotent", {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")

  ws1 <- rMVPA:::as_worker_spec(mspec)
  ws2 <- rMVPA:::as_worker_spec(ws1)

  expect_null(ws2$dataset)
  expect_equal(class(ws2), class(mspec))
})

test_that("strip_dataset and as_worker_spec produce equivalent results", {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")

  stripped <- strip_dataset(mspec)
  worker <- rMVPA:::as_worker_spec(mspec)

  expect_equal(sort(names(stripped)), sort(names(worker)))
  expect_null(stripped$dataset)
  expect_null(worker$dataset)
  expect_equal(class(stripped), class(worker))
})
