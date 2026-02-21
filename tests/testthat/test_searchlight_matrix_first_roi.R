testthat::skip_if_not_installed("neuroim2")

build_matrix_first_mvpa_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 2, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
}

build_matrix_first_rsa_spec <- function() {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 3, blocks = 4)
  dmat <- dist(matrix(rnorm(24 * 12), nrow = 24, ncol = 12))
  rdes <- rsa_design(~ dmat, list(dmat = dmat), block_var = ds$design$block_var)
  rsa_model(ds$dataset, rdes, regtype = "pearson")
}

build_matrix_first_naive_spec <- function() {
  ds <- gen_sample_dataset(
    D = c(4, 4, 4),
    nobs = 24,
    nlevels = 3,
    blocks = 4,
    external_test = TRUE,
    ntest_obs = 24
  )
  naive_xdec_model(ds$dataset, ds$design)
}

run_searchlight_matrix_first <- function(mspec, radius, enabled) {
  force(enabled)
  run_searchlight(mspec, radius = radius, method = "standard", backend = "default")
}

test_that("process_roi_default uses matrix-first input", {
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20, nlevels = 2, blocks = 4)
  vox <- sample(which(ds$dataset$mask > 0), 16)
  roi <- as_roi(data_sample(ds$dataset, vox), ds$dataset)

  mod_spec <- structure(
    list(cv_spec = NULL, crossval = NULL),
    class = c("dummy_model_spec", "model_spec", "list")
  )

  seen_input <- NULL
  fake_train <- function(obj, roi_x, y, indices, ...) {
    seen_input <<- roi_x
    list(ok = TRUE)
  }
  fake_merge <- function(obj, result_set, indices, id, ...) {
    tibble::tibble(
      result = list(NULL),
      indices = list(indices),
      performance = list(NULL),
      id = id,
      error = FALSE,
      error_message = "~"
    )
  }

  old_off <- options(rMVPA.searchlight_mode = "legacy", rMVPA.matrix_first_roi = FALSE)
  on.exit(options(old_off), add = TRUE)
  testthat::with_mocked_bindings(
    rMVPA:::process_roi_default(mod_spec, roi, rnum = 1),
    train_model = fake_train,
    merge_results = fake_merge,
    .package = "rMVPA"
  )
  expect_true(is.matrix(seen_input))
})

test_that("matrix-first ROI path preserves mvpa_model searchlight outputs", {
  mspec <- build_matrix_first_mvpa_spec()

  set.seed(3301)
  res_base <- run_searchlight_matrix_first(mspec, radius = 2, enabled = FALSE)
  set.seed(3301)
  res_fast <- run_searchlight_matrix_first(mspec, radius = 2, enabled = TRUE)

  expect_searchlight_parity(
    reference = res_base,
    candidate = res_fast,
    atol = 1e-10,
    rtol = 1e-8
  )
})

test_that("matrix-first ROI path preserves rsa_model searchlight outputs", {
  mspec <- build_matrix_first_rsa_spec()

  set.seed(3302)
  res_base <- run_searchlight_matrix_first(mspec, radius = 2, enabled = FALSE)
  set.seed(3302)
  res_fast <- run_searchlight_matrix_first(mspec, radius = 2, enabled = TRUE)

  expect_searchlight_parity(
    reference = res_base,
    candidate = res_fast,
    atol = 1e-10,
    rtol = 1e-8
  )
})

test_that("matrix-first ROI path preserves naive_xdec_model searchlight outputs", {
  mspec <- build_matrix_first_naive_spec()

  set.seed(3303)
  res_base <- run_searchlight_matrix_first(mspec, radius = 2, enabled = FALSE)
  set.seed(3303)
  res_fast <- run_searchlight_matrix_first(mspec, radius = 2, enabled = TRUE)

  expect_searchlight_parity(
    reference = res_base,
    candidate = res_fast,
    atol = 1e-10,
    rtol = 1e-8
  )
})
