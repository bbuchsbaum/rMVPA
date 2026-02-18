library(testthat)
library(rMVPA)

futile.logger::flog.threshold(futile.logger::ERROR)

# ---------------------------------------------------------------------------
# Group 1: output_schema dispatch basics
# ---------------------------------------------------------------------------

test_that("output_schema.mvpa_model returns scalar schema for standard model", {
  ds   <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  mdl  <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
  schema <- output_schema(mspec)
  expect_true(is.list(schema))
  expect_equal(names(schema), c("Accuracy", "AUC"))
  expect_true(all(unlist(schema) == "scalar"))
})

test_that("output_schema.default returns NULL for arbitrary list object", {
  obj <- structure(list(x = 1), class = c("some_random_class", "list"))
  expect_null(output_schema(obj))
})

# ---------------------------------------------------------------------------
# Group 2: output_schema.manova_model returns correct schema
# ---------------------------------------------------------------------------

test_that("output_schema.manova_model returns a named list", {
  ds  <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  des <- manova_design(~ block_var + Y, ds$design$train_design)
  mspec <- manova_model(ds$dataset, des)

  schema <- output_schema(mspec)
  expect_true(is.list(schema))
  expect_true(length(schema) > 0)
})

test_that("output_schema.manova_model values are all 'scalar'", {
  ds  <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  des <- manova_design(~ block_var + Y, ds$design$train_design)
  mspec <- manova_model(ds$dataset, des)

  schema <- output_schema(mspec)
  expect_true(all(unlist(schema) == "scalar"))
})

test_that("output_schema.manova_model names match sanitized formula terms", {
  ds  <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  des <- manova_design(~ block_var + Y, ds$design$train_design)
  mspec <- manova_model(ds$dataset, des)

  schema        <- output_schema(mspec)
  expected_names <- rMVPA:::sanitize(labels(terms(des$formula)))
  expect_equal(names(schema), expected_names)
})

test_that("output_schema.manova_model handles single-term formula", {
  ds  <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  des <- manova_design(~ block_var, ds$design$train_design)
  mspec <- manova_model(ds$dataset, des)

  schema <- output_schema(mspec)
  expect_equal(length(schema), 1L)
  expect_equal(names(schema), rMVPA:::sanitize(labels(terms(des$formula))))
  expect_equal(schema[[1]], "scalar")
})

# ---------------------------------------------------------------------------
# Group 3: output_schema.rsa_model returns correct schema
# ---------------------------------------------------------------------------

test_that("output_schema.rsa_model returns a named scalar list", {
  ds     <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 4)
  nobs   <- nrow(ds$design$train_design)
  Dmat   <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))

  rdes  <- rsa_design(~ Dmat, list(Dmat = Dmat))
  rmod  <- rsa_model(ds$dataset, rdes, distmethod = "spearman", regtype = "lm")

  schema <- output_schema(rmod)
  expect_true(is.list(schema))
  expect_equal(length(schema), 1L)
  expect_true(all(unlist(schema) == "scalar"))
  expect_equal(names(schema), "Dmat")
})

test_that("output_schema.rsa_model with two predictors has length 2", {
  ds   <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 4)
  nobs <- nrow(ds$design$train_design)
  D1   <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))
  D2   <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))

  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2))
  rmod <- rsa_model(ds$dataset, rdes, distmethod = "spearman", regtype = "lm")

  schema <- output_schema(rmod)
  expect_equal(length(schema), 2L)
  expect_equal(names(schema), c("D1", "D2"))
  expect_true(all(unlist(schema) == "scalar"))
})

test_that("output_schema.rsa_model names match model_mat names", {
  ds   <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 4)
  nobs <- nrow(ds$design$train_design)
  predictor1 <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))
  predictor2 <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))

  rdes <- rsa_design(~ predictor1 + predictor2,
                     list(predictor1 = predictor1, predictor2 = predictor2))
  rmod <- rsa_model(ds$dataset, rdes, regtype = "lm")

  schema <- output_schema(rmod)
  expect_equal(names(schema), names(rmod$design$model_mat))
})

# ---------------------------------------------------------------------------
# Group 4: build_maps_from_schema unit tests
# ---------------------------------------------------------------------------

test_that("build_maps_from_schema builds correct number of maps for scalar schema", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 20, nlevels = 2)
  schema <- list(accuracy = "scalar", AUC = "scalar")

  ids      <- sample(which(ds$dataset$mask > 0), 10)
  perf_mat <- matrix(runif(20), nrow = 10, ncol = 2,
                     dimnames = list(NULL, c("accuracy", "AUC")))

  maps <- rMVPA:::build_maps_from_schema(schema, perf_mat, ds$dataset, ids)

  expect_equal(length(maps), 2L)
  expect_equal(names(maps), c("accuracy", "AUC"))
})

test_that("build_maps_from_schema scalar maps are spatial objects", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 20, nlevels = 2)
  schema <- list(metric1 = "scalar")

  ids      <- sample(which(ds$dataset$mask > 0), 8)
  perf_mat <- matrix(runif(8), nrow = 8, ncol = 1)

  maps <- rMVPA:::build_maps_from_schema(schema, perf_mat, ds$dataset, ids)

  expect_equal(length(maps), 1L)
  expect_true(!is.null(maps[["metric1"]]))
})

test_that("build_maps_from_schema handles vector schema entries", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 20, nlevels = 2)
  schema <- list(metric1 = "scalar", coefs = "vector[3]")

  ids      <- sample(which(ds$dataset$mask > 0), 10)
  perf_mat <- matrix(runif(40), nrow = 10, ncol = 4)

  maps <- rMVPA:::build_maps_from_schema(schema, perf_mat, ds$dataset, ids)

  expect_equal(length(maps), 4L)
  expect_equal(names(maps), c("metric1", "coefs.1", "coefs.2", "coefs.3"))
})

test_that("build_maps_from_schema vector[1] produces single sub-map", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 20, nlevels = 2)
  schema <- list(v = "vector[1]")

  ids      <- sample(which(ds$dataset$mask > 0), 5)
  perf_mat <- matrix(runif(5), nrow = 5, ncol = 1)

  maps <- rMVPA:::build_maps_from_schema(schema, perf_mat, ds$dataset, ids)

  expect_equal(length(maps), 1L)
  expect_equal(names(maps), "v.1")
})

test_that("build_maps_from_schema errors on unknown schema type", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 20, nlevels = 2)
  schema <- list(bad_metric = "unknown_type")

  ids      <- sample(which(ds$dataset$mask > 0), 10)
  perf_mat <- matrix(runif(10), nrow = 10, ncol = 1)

  expect_error(
    rMVPA:::build_maps_from_schema(schema, perf_mat, ds$dataset, ids),
    "Unknown schema type"
  )
})

test_that("build_maps_from_schema errors on mixed bad type among valid types", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 20, nlevels = 2)
  schema <- list(good = "scalar", bad = "matrix[2x2]")

  ids      <- sample(which(ds$dataset$mask > 0), 5)
  perf_mat <- matrix(runif(10), nrow = 5, ncol = 2)

  expect_error(
    rMVPA:::build_maps_from_schema(schema, perf_mat, ds$dataset, ids),
    "Unknown schema type"
  )
})

# ---------------------------------------------------------------------------
# Group 5: combine_schema_standard activation checks
# ---------------------------------------------------------------------------

test_that("combine_schema_standard is active for standard mvpa_model", {
  ds    <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  mdl   <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")

  expect_false(is.null(output_schema(mspec)))
})

test_that("output_schema is non-NULL for manova_model (schema path active)", {
  ds  <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 2)
  des <- manova_design(~ block_var + Y, ds$design$train_design)
  mspec <- manova_model(ds$dataset, des)

  expect_false(is.null(output_schema(mspec)))
})

test_that("output_schema is non-NULL for rsa_model (schema path active)", {
  ds   <- gen_sample_dataset(c(4, 4, 4), 20, nlevels = 4)
  nobs <- nrow(ds$design$train_design)
  Dmat <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))

  rdes <- rsa_design(~ Dmat, list(Dmat = Dmat))
  rmod <- rsa_model(ds$dataset, rdes, regtype = "lm")

  expect_false(is.null(output_schema(rmod)))
})

# ---------------------------------------------------------------------------
# Group 6: Integration — schema path produces valid searchlight results
# ---------------------------------------------------------------------------

test_that("manova searchlight produces valid results via schema path", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 50, blocks = 3)
  des <- manova_design(~ block_var + Y, ds$design$train_design)
  mspec <- manova_model(ds$dataset, des)

  schema <- output_schema(mspec)
  expect_false(is.null(schema))

  res <- run_searchlight(mspec, radius = 4, method = "standard")

  expect_s3_class(res, "searchlight_result")
  expect_true(length(res$results) > 0)
  expect_equal(names(res$results), names(schema))
})

test_that("manova searchlight result map names match schema names exactly", {
  ds  <- gen_sample_dataset(c(5, 5, 5), 50, blocks = 3)
  des <- manova_design(~ block_var, ds$design$train_design)
  mspec <- manova_model(ds$dataset, des)

  schema <- output_schema(mspec)
  res    <- run_searchlight(mspec, radius = 4, method = "standard")

  expect_equal(sort(names(res$results)), sort(names(schema)))
})

test_that("rsa searchlight produces valid results via schema path", {
  ds   <- gen_sample_dataset(c(5, 5, 5), 50, blocks = 3)
  nobs <- nrow(ds$design$train_design)
  D1   <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))
  D2   <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))

  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2,
                                      block = ds$design$block_var),
                     block_var = "block")
  rmod <- rsa_model(ds$dataset, rdes, regtype = "lm")

  schema <- output_schema(rmod)
  expect_false(is.null(schema))
  expect_equal(names(schema), c("D1", "D2"))

  res <- run_searchlight(rmod, radius = 4, method = "standard")

  expect_s3_class(res, "searchlight_result")
  expect_true(length(res$results) > 0)
  expect_equal(names(res$results), names(schema))
})

test_that("rsa searchlight single predictor produces one result map", {
  ds   <- gen_sample_dataset(c(5, 5, 5), 50, blocks = 3)
  nobs <- nrow(ds$design$train_design)
  Dmat <- dist(matrix(rnorm(nobs * nobs), nobs, nobs))

  rdes <- rsa_design(~ Dmat, list(Dmat = Dmat, block = ds$design$block_var),
                     block_var = "block")
  rmod <- rsa_model(ds$dataset, rdes, regtype = "pearson")

  schema <- output_schema(rmod)
  expect_equal(length(schema), 1L)

  res <- run_searchlight(rmod, radius = 4, method = "standard")

  expect_s3_class(res, "searchlight_result")
  expect_equal(length(res$results), 1L)
  expect_equal(names(res$results), "Dmat")
})

test_that("standard mvpa_model searchlight works via schema path", {
  ds   <- gen_sample_dataset(c(5, 5, 5), 30, nlevels = 2, blocks = 3)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl  <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = cval)

  schema <- output_schema(mspec)
  expect_false(is.null(schema))

  res <- run_searchlight(mspec, radius = 4, method = "standard")

  expect_s3_class(res, "searchlight_result")
  expect_equal(sort(names(res$results)), sort(names(schema)))
})
