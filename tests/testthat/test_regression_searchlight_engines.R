# Regression tests for searchlight engine selection, argument validation,
# non-finite data handling, error attribution and the geometry cache.
# Each test here failed before the fixes and passes after them.
testthat::skip_if_not_installed("neuroim2")

reg_gen <- function(seed, nlevels = 2, nobs = 32, na_cols = 0) {
  set.seed(seed)
  gen_sample_dataset(c(5, 5, 5), nobs = nobs, nlevels = nlevels, blocks = 4,
                     na_cols = na_cols)
}

reg_dual_lda_spec <- function(seed = 2101, na_cols = 0, feature_selector = NULL) {
  ds <- reg_gen(seed, na_cols = na_cols)
  mvpa_model(load_model("dual_lda"), ds$dataset, ds$design, "classification",
             crossval = blocked_cross_validation(ds$design$block_var),
             tune_grid = data.frame(gamma = 1e-2),
             feature_selector = feature_selector)
}

reg_corclass_spec <- function(seed = 2102) {
  ds <- reg_gen(seed, nlevels = 3, nobs = 48)
  mvpa_model(load_model("corclass"), ds$dataset, ds$design, "classification",
             crossval = blocked_cross_validation(ds$design$block_var))
}

reg_maps <- function(res) lapply(res$results, function(m) as.numeric(neuroim2::values(m)))

quiet_run <- function(...) suppressWarnings(suppressMessages(run_searchlight(...)))

reg_error_message <- function(expr) {
  tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
}

test_that("dual_lda with a feature_selector under auto declines and matches legacy", {
  mspec <- reg_dual_lda_spec(feature_selector = feature_selector("FTest", "top_k", 10))

  auto <- quiet_run(mspec, radius = 2, engine = "auto", backend = "default")
  legacy <- quiet_run(mspec, radius = 2, engine = "legacy", backend = "default")

  expect_false(identical(attr(auto, "searchlight_engine"), "dual_lda_fast"))
  expect_identical(attr(auto, "searchlight_engine"), "legacy")
  expect_identical(sort(names(auto$results)), sort(names(legacy$results)))
  expect_identical(reg_maps(auto), reg_maps(legacy))
})

test_that("dual_lda with non-finite data under auto falls back to legacy and matches it", {
  mspec <- reg_dual_lda_spec(seed = 2103, na_cols = 3)

  auto <- quiet_run(mspec, radius = 2, engine = "auto", backend = "default")
  legacy <- quiet_run(mspec, radius = 2, engine = "legacy", backend = "default")

  expect_identical(attr(auto, "searchlight_engine"), "legacy")
  expect_identical(sort(names(auto$results)), sort(names(legacy$results)))
  a <- reg_maps(auto)
  b <- reg_maps(legacy)
  for (m in names(b)) {
    expect_identical(is.na(a[[m]]), is.na(b[[m]]), info = m)
    expect_identical(a[[m]], b[[m]], info = m)
  }
})

test_that("explicit dual_lda_fast with non-finite data errors clearly", {
  mspec <- reg_dual_lda_spec(seed = 2104, na_cols = 3)
  msg <- reg_error_message(quiet_run(mspec, radius = 2, engine = "dual_lda_fast",
                                     backend = "default"))
  expect_false(is.na(msg))
  expect_match(msg, "dual_lda_fast", fixed = TRUE)
  expect_match(msg, "non-finite", fixed = TRUE)
})

test_that("radius validation errors identically under auto, legacy and fast engines", {
  mspec <- reg_corclass_spec()
  msg_legacy <- reg_error_message(quiet_run(mspec, radius = 0.5, engine = "legacy",
                                            backend = "default"))
  expect_match(msg_legacy, "outside allowable range", fixed = TRUE)
  for (eng in c("auto", "aggregate_fast")) {
    msg <- reg_error_message(quiet_run(mspec, radius = 0.5, engine = eng,
                                       backend = "default"))
    expect_identical(msg, msg_legacy, info = eng)
  }
})

test_that("string combiner validation errors under auto as under legacy", {
  mspec <- reg_corclass_spec()
  msg_legacy <- reg_error_message(quiet_run(mspec, radius = 2, engine = "legacy",
                                            combiner = "pool", backend = "default"))
  expect_match(msg_legacy, "Unknown string combiner", fixed = TRUE)
  msg_auto <- reg_error_message(quiet_run(mspec, radius = 2, engine = "auto",
                                          combiner = "pool", backend = "default"))
  expect_identical(msg_auto, msg_legacy)
})

test_that("a custom standard combiner keeps the fast engine off under auto", {
  # Standard-method mvpa_model results are assembled by the output schema,
  # so the combiner itself is not called here on either path; the check is
  # that auto does not hand the run to a fast engine that ignores it.
  mspec <- reg_corclass_spec()
  counting_combiner <- function(model_spec, good_results, bad_results) {
    rMVPA:::combine_standard(model_spec, good_results, bad_results)
  }

  auto <- quiet_run(mspec, radius = 2, engine = "auto", combiner = counting_combiner,
                    backend = "default")
  expect_identical(attr(auto, "searchlight_engine"), "legacy")

  legacy <- quiet_run(mspec, radius = 2, engine = "legacy", combiner = counting_combiner,
                      backend = "default")
  expect_identical(reg_maps(auto), reg_maps(legacy))
})

test_that("a custom randomized combiner is called under auto and matches legacy", {
  # The randomized legacy path calls the combiner; the dual_lda sampled fast
  # path only supports the built-in average. Under auto the run must reach
  # the general path, so the closure is called and no engine-failure
  # fallback warning is emitted.
  mspec <- reg_dual_lda_spec(seed = 2108)
  calls <- 0L
  counting_combiner <- function(model_spec, good_results, bad_results, ...) {
    calls <<- calls + 1L
    rMVPA:::combine_randomized(model_spec, good_results, bad_results, ...)
  }

  run_rand <- function(engine) {
    set.seed(7)
    msgs <- character(0)
    res <- withCallingHandlers(
      suppressMessages(run_searchlight(mspec, radius = 2, method = "randomized",
                                       niter = 2, engine = engine,
                                       combiner = counting_combiner,
                                       backend = "default")),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    list(res = res, msgs = msgs)
  }

  auto <- run_rand("auto")
  expect_gt(calls, 0L)
  expect_identical(attr(auto$res, "searchlight_engine"), "legacy")
  expect_false(any(grepl("searchlight engine", auto$msgs, fixed = TRUE)))

  calls_before_legacy <- calls
  legacy <- run_rand("legacy")
  expect_gt(calls, calls_before_legacy)
  expect_identical(reg_maps(auto$res), reg_maps(legacy$res))
})

test_that("geometry cache does not return another mask's centres on a signature collision", {
  base_ds <- reg_gen(2105)
  space_obj <- neuroim2::space(base_ds$dataset$mask)
  make_dset <- function(idx) {
    arr <- array(FALSE, dim = dim(space_obj)[1:3])
    arr[idx] <- TRUE
    mask <- neuroim2::LogicalNeuroVol(arr, space_obj)
    mvpa_dataset(train_data = base_ds$dataset$train_data, mask = mask)
  }
  idx1 <- c(10, 20, 30, 40, 50, 60)
  idx2 <- c(11, 18, 31, 40, 50, 60)

  rMVPA:::.searchlight_geometry_cache_clear()
  d1 <- make_dset(idx1)
  d2 <- make_dset(idx2)
  ids <- function(sl) vapply(sl, function(w) as.integer(w@parent_index), integer(1))

  expect_setequal(ids(get_searchlight(d1, type = "standard", radius = 2)), idx1)
  expect_setequal(ids(get_searchlight(d2, type = "standard", radius = 2)), idx2)
})

test_that("ineligibility messages name the engine that declined", {
  set.seed(2106)
  ds <- gen_sample_dataset(c(5, 5, 5), 32, blocks = 4, na_cols = 3)
  D1 <- dist(matrix(rnorm(32 * 4), 32))
  D2 <- dist(matrix(rnorm(32 * 4), 32))
  rdes <- rsa_design(~ D1 + D2, list(D1 = D1, D2 = D2, block = ds$design$block_var),
                     block_var = "block")
  rms <- rsa_model(ds$dataset, rdes, distmethod = "pearson", regtype = "pearson",
                   check_collinearity = FALSE)
  msg <- reg_error_message(quiet_run(rms, radius = 2, engine = "rsa_fast",
                                     backend = "default"))
  expect_match(msg, "rsa_fast: data contain missing or non-finite values", fixed = TRUE)
  expect_no_match(msg, "aggregate_fast", fixed = TRUE)

  mspec <- reg_corclass_spec(2107)
  mspec$dataset$train_data <- neuroim2::NeuroVec(
    {
      arr <- neuroim2::as.array(mspec$dataset$train_data)
      arr[1, 1, 1, ] <- NA
      arr
    },
    neuroim2::space(mspec$dataset$train_data)
  )
  msg_agg <- reg_error_message(quiet_run(mspec, radius = 2, engine = "aggregate_fast",
                                         backend = "default"))
  expect_match(msg_agg, "aggregate_fast: data contain", fixed = TRUE)
})
