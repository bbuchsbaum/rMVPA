testthat::skip_if_not_installed("neuroim2")

# Regression tests for engine selection under engine = "auto": a fast engine may
# be chosen only when it computes the estimator named in the model spec. Before
# this fix, every multiclass mvpa_model was routed to SWIFT, which replaced the
# requested classifier with SWIFT's own nearest-class-mean estimator.

build_identity_mspec <- function(model_name, nlevels = 3L, D = c(4, 4, 4),
                                 nobs = 48, blocks = 4) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = nlevels, blocks = blocks)
  cval <- blocked_cross_validation(ds$design$block_var)
  mvpa_model(
    model = suppressWarnings(load_model(model_name)),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cval
  )
}

identity_map_values <- function(res) {
  lapply(res$results, function(m) as.numeric(neuroim2::values(m)))
}

test_that("auto does not resolve multiclass classifiers to SWIFT", {
  for (model_name in c("corclass", "sda_notune", "svmLinear")) {
    set.seed(5101)
    mspec <- build_identity_mspec(model_name)

    # SWIFT is still registered and eligible for an explicit request ...
    tbl <- searchlight_engines(mspec, method = "standard")
    expect_true(tbl$eligible[tbl$engine == "swift"], info = model_name)

    # ... but auto never selects it.
    expect_identical(
      .resolve_searchlight_engine(mspec, method = "standard", engine = "auto"),
      "legacy",
      info = model_name
    )
    explained <- explain_searchlight_engine(mspec, method = "standard", engine = "auto")
    expect_identical(explained$engine[explained$selected], "legacy", info = model_name)
  }
})

test_that("auto and legacy give identical corclass searchlight maps", {
  set.seed(5102)
  mspec <- build_identity_mspec("corclass")

  set.seed(1)
  res_auto <- run_searchlight(mspec, radius = 2, method = "standard",
                              backend = "default")
  set.seed(1)
  res_legacy <- run_searchlight(mspec, radius = 2, method = "standard",
                                backend = "default", engine = "legacy")

  expect_identical(attr(res_auto, "searchlight_engine"), "legacy")
  expect_false("SWIFT_Info" %in% names(res_auto$results))
  expect_identical(sort(names(res_auto$results)), sort(names(res_legacy$results)))
  expect_equal(identity_map_values(res_auto)[names(res_legacy$results)],
               identity_map_values(res_legacy))
})

test_that("dual_lda keeps its exact fast engine under auto", {
  set.seed(5103)
  mspec <- build_identity_mspec("dual_lda")
  expect_identical(
    .resolve_searchlight_engine(mspec, method = "standard", engine = "auto"),
    "dual_lda_fast"
  )
})

test_that("explicit swift request announces the estimator substitution", {
  set.seed(5104)
  mspec <- build_identity_mspec("corclass")
  expect_message(
    res <- run_searchlight(mspec, radius = 2, method = "standard",
                           backend = "default", engine = "swift"),
    regexp = "not the requested 'corclass' classifier"
  )
  expect_identical(attr(res, "searchlight_engine"), "swift")
  expect_identical(attr(res, "searchlight_estimator"), "swift_nearest_mean")
})
