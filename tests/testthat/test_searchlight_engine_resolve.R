library(rMVPA)

# Regression tests for searchlight engine resolution/execution.
#
# 1. explain_searchlight_engine() must report the engine that would actually
#    run, including for model classes that do NOT inherit "mvpa_model" (e.g.
#    naive_xdec_model). Previously it hardcoded an `inherits(., "mvpa_model")`
#    gate and always reported "legacy" for such models, contradicting both
#    searchlight_engines() and the engine that actually executed.
# 2. .execute_searchlight_engine() must fail loudly on an engine it cannot run,
#    rather than returning NULL (which the caller treats as a successful empty
#    result).

make_eligible_naive_xdec <- function() {
  toy <- gen_sample_dataset(
    D = c(5, 5, 5), nobs = 60, nlevels = 3, blocks = 3, external_test = TRUE
  )
  # return_predictions = FALSE is required for the fast-metric kernel path.
  naive_xdec_model(toy$dataset, toy$design, link_by = NULL,
                   return_predictions = FALSE)
}

test_that("explain_searchlight_engine reports the fast path for naive_xdec_model", {
  ms <- make_eligible_naive_xdec()
  skip_if_not(isTRUE(rMVPA:::.is_naive_xdec_fast_path(ms, "standard")),
              "naive_xdec fast path not eligible in this environment")

  ex <- explain_searchlight_engine(ms, method = "standard", engine = "auto")
  selected <- ex$engine[ex$selected]

  expect_identical(selected, "naive_xdec_fast")
  # The eligible flag and the selected flag must agree with each other.
  expect_true(ex$eligible[ex$engine == "naive_xdec_fast"])
})

test_that("explain_searchlight_engine matches the engine that actually runs", {
  ms <- make_eligible_naive_xdec()
  skip_if_not(isTRUE(rMVPA:::.is_naive_xdec_fast_path(ms, "standard")),
              "naive_xdec fast path not eligible in this environment")

  ex <- explain_searchlight_engine(ms, method = "standard", engine = "auto")
  selected <- ex$engine[ex$selected]

  suppressWarnings(
    res <- run_searchlight(ms, radius = 2, method = "standard")
  )
  expect_identical(attr(res, "searchlight_engine"), selected)
})

test_that("explain_searchlight_engine honours explicit and ineligible requests", {
  ms <- make_eligible_naive_xdec()
  skip_if_not(isTRUE(rMVPA:::.is_naive_xdec_fast_path(ms, "standard")),
              "naive_xdec fast path not eligible in this environment")

  explicit <- explain_searchlight_engine(ms, method = "standard",
                                         engine = "naive_xdec_fast")
  expect_identical(explicit$engine[explicit$selected], "naive_xdec_fast")

  # An explicit legacy request always resolves to legacy.
  legacy_req <- explain_searchlight_engine(ms, method = "standard",
                                           engine = "legacy")
  expect_identical(legacy_req$engine[legacy_req$selected], "legacy")

  # The fast path is standard-only; a randomized run falls back to legacy.
  randomized <- explain_searchlight_engine(ms, method = "randomized",
                                           engine = "auto")
  expect_identical(randomized$engine[randomized$selected], "legacy")
})

test_that(".resolve_searchlight_engine defaults to legacy for unknown classes", {
  spec <- structure(list(), class = c("some_other_model", "model_spec", "list"))
  expect_identical(
    rMVPA:::.resolve_searchlight_engine(spec, method = "standard", engine = "auto"),
    "legacy"
  )
})

test_that(".execute_searchlight_engine stops on an unhandled engine", {
  # naive_xdec_fast is registered but has no executor branch in
  # .execute_searchlight_engine (it dispatches via its own S3 method instead).
  # Reaching this executor with such an engine is a wiring bug and must error
  # rather than silently return NULL.
  expect_error(
    rMVPA:::.execute_searchlight_engine(
      engine = "naive_xdec_fast",
      model_spec = structure(list(), class = "mvpa_model"),
      radius = 2,
      method = "standard"
    ),
    "no executor for searchlight engine"
  )
})
