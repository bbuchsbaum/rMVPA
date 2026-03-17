test_that("load_model checks for required packages", {
  # Register a model that requires a non-existent package
  fake_spec <- list(
    type = "Classification",
    library = "this_package_does_not_exist_12345",
    label = "fake_model",
    parameters = data.frame(parameter = "x", class = "numeric", label = "X"),
    grid = function(x, y, len = NULL) data.frame(x = 1),
    fit = function(x, y, ...) NULL,
    predict = function(modelFit, newdata, ...) NULL,
    prob = function(modelFit, newdata, ...) NULL
  )
  register_mvpa_model("fake_model", fake_spec)
  on.exit(rm("fake_model", envir = MVPAModels), add = TRUE)

  expect_warning(
    load_model("fake_model"),
    "requires package.*not installed"
  )
  expect_warning(
    load_model("fake_model"),
    "install.packages"
  )
})


test_that("load_model succeeds for models with available packages", {
  # corclass has no external library requirement
  mod <- load_model("corclass")
  expect_true(is.list(mod))
  expect_equal(mod$label, "corclass")
})


test_that("require_package gives helpful error message", {
  expect_error(
    rMVPA:::require_package("this_package_does_not_exist_12345",
                            "for testing purposes"),
    "this_package_does_not_exist_12345.*required.*for testing purposes"
  )
  expect_error(
    rMVPA:::require_package("this_package_does_not_exist_12345"),
    "install.packages"
  )
})


test_that("require_package succeeds for installed packages", {
  # stats is always available
  expect_silent(rMVPA:::require_package("stats"))
})


test_that("glmnet_opt fit uses namespace-qualified glmnet call", {
  fit_body <- paste(deparse(body(MVPAModels$glmnet_opt$fit)), collapse = "\n")
  expect_match(fit_body, "glmnet::glmnet\\(")
})
