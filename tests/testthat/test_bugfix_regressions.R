# Regression tests for bug fixes
# These tests ensure that previously fixed bugs don't regress

context("bugfix regressions")

# =============================================================================
# Test: predicted_class() returns a value (not NULL)
# Bug: The function assigned result to 'pclass' but never returned it
# Fixed in: R/performance.R
# =============================================================================
test_that("predicted_class returns the predicted class labels", {
  # Create a simple probability matrix
  prob <- matrix(c(0.2, 0.8,
                   0.6, 0.4,
                   0.1, 0.9),
                 nrow = 3, byrow = TRUE,
                 dimnames = list(NULL, c("A", "B")))

  result <- predicted_class(prob)

  # Should return a character vector, not NULL
  expect_false(is.null(result))
  expect_true(is.character(result))
  expect_equal(length(result), 3)

  # Check the predictions are correct (highest probability class)
  expect_equal(result, c("B", "A", "B"))
})

test_that("predicted_class handles ties correctly", {
  # Create matrix with equal probabilities
  prob <- matrix(c(0.5, 0.5,
                   0.3, 0.7),
                 nrow = 2, byrow = TRUE,
                 dimnames = list(NULL, c("X", "Y")))

  result <- predicted_class(prob)

  # Should still return valid class labels
  expect_false(is.null(result))
  expect_true(all(result %in% c("X", "Y")))
})

# =============================================================================
# Test: regional_mvpa_result() works with default fits=NULL argument
# Bug: Default argument was fits=fits (self-referential)
# Fixed in: R/regional.R
# =============================================================================
test_that("regional_mvpa_result works with default fits argument", {
  # Create minimal mock objects
  model_spec <- list(dataset = "mock_dataset")
  performance_table <- data.frame(accuracy = c(0.8, 0.85))
  prediction_table <- data.frame(
    observed = factor(c("A", "B", "A", "B")),
    predicted = factor(c("A", "B", "B", "B"))
  )
  vol_results <- list(accuracy = "mock_vol")

  # This should work without specifying fits (defaults to NULL)
  result <- regional_mvpa_result(
    model_spec = model_spec,
    performance_table = performance_table,
    prediction_table = prediction_table,
    vol_results = vol_results
  )

  expect_true(inherits(result, "regional_mvpa_result"))
  expect_null(result$fits)
  expect_equal(result$model_spec, model_spec)
})

test_that("regional_mvpa_result works with explicit fits argument", {
  model_spec <- list(dataset = "mock_dataset")
  performance_table <- data.frame(accuracy = 0.9)
  prediction_table <- data.frame(observed = factor("A"), predicted = factor("A"))
  vol_results <- list()
  fits <- list(fit1 = "mock_fit")

  result <- regional_mvpa_result(
    model_spec = model_spec,
    performance_table = performance_table,
    prediction_table = prediction_table,
    vol_results = vol_results,
    fits = fits
  )

  expect_equal(result$fits, fits)
})

# =============================================================================
# Test: combine_prediction_tables() returns a value (not NULL)
# Bug: Function created prediction_table but didn't return it
# Fixed in: R/regional.R
# =============================================================================
test_that("combine_prediction_tables returns a prediction table", {
  # Create sample prediction tables
  observed <- factor(c("A", "B", "A", "B"))

  predtab1 <- data.frame(
    .rownum = 1:4,
    roinum = rep(1, 4),
    observed = observed,
    prob_A = c(0.8, 0.3, 0.7, 0.2),
    prob_B = c(0.2, 0.7, 0.3, 0.8)
  )

  predtab2 <- data.frame(
    .rownum = 1:4,
    roinum = rep(2, 4),
    observed = observed,
    prob_A = c(0.6, 0.4, 0.9, 0.1),
    prob_B = c(0.4, 0.6, 0.1, 0.9)
  )

  result <- combine_prediction_tables(list(predtab1, predtab2))

  # Should return a tibble/data.frame, not NULL
  expect_false(is.null(result))
  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true(".rownum" %in% names(result))
  expect_true("observed" %in% names(result))
  expect_true("predicted" %in% names(result))
})

test_that("combine_prediction_tables handles collapse_regions=TRUE", {
  observed <- factor(c("A", "B"))

  predtab1 <- data.frame(
    .rownum = 1:2,
    roinum = rep(1, 2),
    observed = observed,
    prob_A = c(0.8, 0.3),
    prob_B = c(0.2, 0.7)
  )

  predtab2 <- data.frame(
    .rownum = 1:2,
    roinum = rep(2, 2),
    observed = observed,
    prob_A = c(0.6, 0.4),
    prob_B = c(0.4, 0.6)
  )

  result <- combine_prediction_tables(list(predtab1, predtab2), collapse_regions = TRUE)

  expect_false(is.null(result))
  expect_true(is.data.frame(result))
  # With collapse_regions=TRUE, should have fewer rows (grouped by .rownum)
  expect_equal(nrow(result), 2)
})

test_that("combine_prediction_tables throws error for unsupported types", {
  # Test with numeric observed (should error - not implemented)
  predtab <- data.frame(
    .rownum = 1:3,
    roinum = rep(1, 3),
    observed = c(1.0, 2.0, 3.0),  # numeric, not factor
    prob_A = c(0.5, 0.5, 0.5)
  )

  expect_error(
    combine_prediction_tables(list(predtab)),
    "combining continuous predictions not implemented"
  )
})

test_that("combine_prediction_tables throws error for invalid types",
{
  # Test with logical observed (neither factor/character nor numeric)
  predtab <- data.frame(
    .rownum = 1:2,
    roinum = rep(1, 2),
    observed = c(TRUE, FALSE),
    prob_A = c(0.5, 0.5)
  )

  expect_error(
    combine_prediction_tables(list(predtab)),
    "observed values must be character, factor, or numeric"
  )
})

# =============================================================================
# Test: merge_results.regional_mvpa_result() extracts prediction tables
# Bug: Function passed regional_mvpa_result objects instead of prediction tables
# Fixed in: R/regional.R
# =============================================================================
test_that("merge_results.regional_mvpa_result extracts prediction tables correctly", {
  # Create two regional_mvpa_result objects with prediction tables
  observed <- factor(c("A", "B", "A", "B"))

  predtab1 <- data.frame(
    .rownum = 1:4,
    roinum = rep(1, 4),
    observed = observed,
    prob_A = c(0.8, 0.3, 0.7, 0.2),
    prob_B = c(0.2, 0.7, 0.3, 0.8)
  )

  predtab2 <- data.frame(
    .rownum = 1:4,
    roinum = rep(2, 4),
    observed = observed,
    prob_A = c(0.6, 0.4, 0.9, 0.1),
    prob_B = c(0.4, 0.6, 0.1, 0.9)
  )

  result1 <- regional_mvpa_result(
    model_spec = list(),
    performance_table = data.frame(),
    prediction_table = predtab1,
    vol_results = list()
  )

  result2 <- regional_mvpa_result(
    model_spec = list(),
    performance_table = data.frame(),
    prediction_table = predtab2,
    vol_results = list()
  )

  # This should NOT error - previously it passed the wrong type
  merged <- merge_results(result1, result2)

  expect_false(is.null(merged))
  expect_true(is.data.frame(merged))
  expect_true(nrow(merged) > 0)
})

test_that("merge_results.regional_mvpa_result handles NULL prediction tables", {
  # Create results with NULL prediction tables
  result1 <- regional_mvpa_result(
    model_spec = list(),
    performance_table = data.frame(),
    prediction_table = NULL,
    vol_results = list()
  )

  result2 <- regional_mvpa_result(
    model_spec = list(),
    performance_table = data.frame(),
    prediction_table = NULL,
    vol_results = list()
  )

  # Should warn and return empty tibble
  expect_warning(
    merged <- merge_results(result1, result2),
    "No valid prediction tables to merge"
  )

  expect_true(is.data.frame(merged))
  expect_equal(nrow(merged), 0)
})

# =============================================================================
# Test: predict.list_model() returns a value (not NULL)
# Bug: Function assigned result to 'res' but never returned it
# Fixed in: R/model_fit.R
# =============================================================================
test_that("predict.list_model returns predictions", {
  # Create a simple mock list_model
  # We'll use basic linear models for simplicity
  set.seed(123)
  train_data <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
  y <- train_data$x1 + train_data$x2 + rnorm(20, sd = 0.1)

  # Fit two simple lm models
  fit1 <- lm(y ~ x1, data = train_data)
  fit2 <- lm(y ~ x2, data = train_data)

  # Create list_model
  list_mod <- list(fit1, fit2)
  class(list_mod) <- c("list_model", "list")

  # Create new data for prediction
  newdata <- data.frame(x1 = c(0.5, -0.5), x2 = c(0.5, -0.5))

  result <- predict(list_mod, newdata)

  # Should return a list of predictions, not NULL
  expect_false(is.null(result))
  expect_true(is.list(result))
  expect_equal(length(result), 2)

  # Each element should have predictions
  expect_equal(length(result[[1]]), 2)
  expect_equal(length(result[[2]]), 2)
})

test_that("predict.list_model errors on NULL newdata", {
  list_mod <- list()
  class(list_mod) <- c("list_model", "list")

  expect_error(
    predict(list_mod, newdata = NULL),
    "newdata cannot be null"
  )
})

# =============================================================================
# Integration test: Full regional analysis pipeline
# This ensures all the fixes work together
# =============================================================================
test_that("full regional analysis pipeline works with fixes", {
  skip_if_not_installed("neuroim2")

  # Generate sample dataset
  dset <- gen_sample_dataset(
    c(6, 6, 4),
    nobs = 50,
    nlevels = 2,
    data_mode = "image",
    response_type = "categorical"
  )

  # Create cross-validation
  cval <- blocked_cross_validation(dset$design$block_var)

  # Create region mask with 3 ROIs
  region_mask <- neuroim2::NeuroVol(
    sample(1:3, size = length(dset$dataset$mask), replace = TRUE),
    space = neuroim2::space(dset$dataset$mask)
  )

  # Load model and create model spec
  model <- load_model("sda_notune")
  mspec <- mvpa_model(
    model,
    dset$dataset,
    dset$design,
    model_type = "classification",
    crossval = cval,
    return_predictions = TRUE
  )

  # Run regional analysis
  res <- run_regional(mspec, region_mask)

  # Verify result structure
  expect_true(inherits(res, "regional_mvpa_result"))
  expect_false(is.null(res$prediction_table))
  expect_true(nrow(res$prediction_table) > 0)

  # Verify predicted_class is used correctly in the pipeline
  if (!is.null(res$prediction_table$predicted)) {
    expect_true(is.factor(res$prediction_table$predicted) ||
                is.character(res$prediction_table$predicted))
  }
})
