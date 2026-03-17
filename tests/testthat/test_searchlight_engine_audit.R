test_that("run_searchlight.default stamps engine attr for fast-path execution", {
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
        engine = "swift"
      )
    },
    .package = "rMVPA"
  )

  expect_identical(attr(out, "searchlight_engine"), "swift")
})

test_that("run_searchlight.default stamps legacy engine attr for fallback", {
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
        engine = "legacy"
      )
    },
    run_searchlight_base = function(...) {
      structure(list(results = list()), class = "searchlight_result")
    },
    .package = "rMVPA"
  )

  expect_identical(attr(out, "searchlight_engine"), "legacy")
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

test_that("auto engine failure falls back to legacy with warning", {
  warns <- character(0)
  ret <- withCallingHandlers(
    testthat::with_mocked_bindings(
      rMVPA:::`.run_searchlight_engine.mvpa_model`(
        model_spec = structure(list(), class = "mvpa_model"),
        radius = 2,
        method = "standard",
        engine = "auto"
      ),
      `.resolve_searchlight_engine.mvpa_model` = function(...) "swift",
      .execute_searchlight_engine = function(...) stop("synthetic engine failure"),
      .package = "rMVPA"
    ),
    warning = function(w) {
      warns[[length(warns) + 1L]] <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )
  expect_false(ret$handled)
  expect_identical(ret$engine, "legacy")
  expect_true(any(grepl("falling back to legacy", warns)))
})
