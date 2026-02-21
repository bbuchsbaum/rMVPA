test_that("run_searchlight.default records engine audit metadata for fast-path execution", {
  out <- testthat::with_mocked_bindings(
    rMVPA:::run_searchlight.default(
      model_spec = list(),
      radius = 2,
      method = "standard"
    ),
    .run_searchlight_engine = function(...) {
      list(
        handled = TRUE,
        result = structure(
          list(results = list()),
          class = "searchlight_result"
        ),
        engine = "swift",
        requested = "auto",
        reason = "auto_swift"
      )
    },
    .package = "rMVPA"
  )

  expect_identical(attr(out, "searchlight_engine"), "swift")
  expect_identical(attr(out, "searchlight_engine_requested"), "auto")
  expect_identical(attr(out, "searchlight_engine_resolved"), "swift")
  expect_identical(attr(out, "searchlight_engine_reason"), "auto_swift")
})

test_that("run_searchlight.default records engine audit metadata for legacy fallback", {
  out <- testthat::with_mocked_bindings(
    rMVPA:::run_searchlight.default(
      model_spec = list(),
      radius = 2,
      method = "standard"
    ),
    .run_searchlight_engine = function(...) {
      list(
        handled = FALSE,
        result = NULL,
        engine = "legacy",
        requested = "auto",
        reason = "auto_ineligible"
      )
    },
    run_searchlight_base = function(...) {
      structure(list(results = list()), class = "searchlight_result")
    },
    .package = "rMVPA"
  )

  expect_identical(attr(out, "searchlight_engine"), "legacy")
  expect_identical(attr(out, "searchlight_engine_requested"), "auto")
  expect_identical(attr(out, "searchlight_engine_resolved"), "legacy")
  expect_identical(attr(out, "searchlight_engine_reason"), "auto_ineligible")
})

test_that("run_searchlight.default records fallback error context", {
  out <- testthat::with_mocked_bindings(
    rMVPA:::run_searchlight.default(
      model_spec = list(),
      radius = 2,
      method = "standard"
    ),
    .run_searchlight_engine = function(...) {
      list(
        handled = FALSE,
        result = NULL,
        engine = "legacy",
        requested = "auto",
        reason = "fast_failed_fallback",
        fallback_error = "synthetic failure"
      )
    },
    run_searchlight_base = function(...) {
      structure(list(results = list()), class = "searchlight_result")
    },
    .package = "rMVPA"
  )

  expect_identical(attr(out, "searchlight_engine_reason"), "fast_failed_fallback")
  expect_identical(attr(out, "searchlight_engine_fallback_error"), "synthetic failure")
})

test_that("strict engine request + execution failure stops with error", {
  expect_error(
    testthat::with_mocked_bindings(
      rMVPA:::`.run_searchlight_engine.mvpa_model`(
        model_spec = structure(list(), class = "mvpa_model"),
        radius = 2,
        method = "standard",
        engine = "swift"
      ),
      `.resolve_searchlight_engine.mvpa_model` = function(...) "swift",
      .execute_searchlight_engine = function(...) stop("synthetic engine failure"),
      .package = "rMVPA"
    ),
    "Requested searchlight engine 'swift' failed"
  )
})

test_that("fail_open = FALSE + auto engine failure stops with error", {
  expect_error(
    testthat::with_mocked_bindings(
      rMVPA:::`.run_searchlight_engine.mvpa_model`(
        model_spec = structure(list(), class = "mvpa_model"),
        radius = 2,
        method = "standard",
        engine = "auto"
      ),
      `.resolve_searchlight_engine.mvpa_model` = function(...) "swift",
      .execute_searchlight_engine = function(...) stop("synthetic engine failure"),
      .searchlight_fail_open_enabled = function() FALSE,
      .package = "rMVPA"
    ),
    "fail-open is disabled"
  )
})
