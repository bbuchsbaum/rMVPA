with_backend_policy_options <- function(code) {
  keys <- c(
    "rMVPA.searchlight_backend_default",
    "rMVPA.searchlight_profile",
    "rMVPA.searchlight_mode"
  )
  old <- options()[keys]
  on.exit(options(old), add = TRUE)
  force(code)
}

test_that("searchlight backend resolver defaults to auto and ignores ambient options", {
  with_backend_policy_options({
    options(
      rMVPA.searchlight_backend_default = "default",
      rMVPA.searchlight_profile = "legacy",
      rMVPA.searchlight_mode = "legacy"
    )

    got <- rMVPA:::.resolve_searchlight_backend(
      backend = c("default", "shard", "auto"),
      missing_backend = TRUE
    )
    expect_identical(got, "auto")
  })
})

test_that("searchlight backend resolver keeps explicit backend override", {
  with_backend_policy_options({
    got <- rMVPA:::.resolve_searchlight_backend(
      backend = "default",
      missing_backend = FALSE
    )
    expect_identical(got, "default")
  })
})

test_that("run_searchlight.default applies backend policy only when backend is omitted", {
  with_backend_policy_options({
    fake_base <- function(model_spec, radius, method, niter, combiner, drop_probs,
                          fail_fast, k, backend, verbose, ...) {
      backend
    }

    out_missing <- testthat::with_mocked_bindings(
      rMVPA:::run_searchlight.default(
        model_spec = list(),
        radius = 2,
        method = "standard"
      ),
      run_searchlight_base = fake_base,
      .package = "rMVPA"
    )
    expect_identical(out_missing, "auto")

    out_explicit <- testthat::with_mocked_bindings(
      rMVPA:::run_searchlight.default(
        model_spec = list(),
        radius = 2,
        method = "standard",
        backend = "default"
      ),
      run_searchlight_base = fake_base,
      .package = "rMVPA"
    )
    expect_identical(out_explicit, "default")
  })
})
