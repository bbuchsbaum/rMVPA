with_mode_options <- function(code) {
  keys <- c(
    "rMVPA.searchlight_mode",
    "rMVPA.searchlight_profile",
    "rMVPA.fold_cache_enabled",
    "rMVPA.searchlight_geometry_cache",
    "rMVPA.clustered_nn_fastpath",
    "rMVPA.matrix_first_roi",
    "rMVPA.fast_filter_roi",
    "rMVPA.rsa_fast_kernel",
    "rMVPA.naive_xdec_fast_kernel",
    "rMVPA.searchlight_backend_default"
  )
  old <- options()[keys]
  on.exit(options(old), add = TRUE)
  force(code)
}

test_that("searchlight mode is deterministic fast", {
  with_mode_options({
    options(rMVPA.searchlight_mode = "legacy", rMVPA.searchlight_profile = "legacy")
    expect_identical(rMVPA:::.resolve_searchlight_mode(), "fast")
    expect_identical(rMVPA::searchlight_mode(), "fast")

    rMVPA::searchlight_mode("legacy")
    expect_identical(rMVPA::searchlight_mode(), "fast")
  })
})

test_that("stable fast paths are always enabled", {
  with_mode_options({
    options(
      rMVPA.searchlight_mode = "legacy",
      rMVPA.searchlight_profile = "legacy",
      rMVPA.fold_cache_enabled = FALSE,
      rMVPA.searchlight_geometry_cache = FALSE,
      rMVPA.clustered_nn_fastpath = FALSE,
      rMVPA.matrix_first_roi = FALSE,
      rMVPA.fast_filter_roi = FALSE,
      rMVPA.rsa_fast_kernel = FALSE,
      rMVPA.naive_xdec_fast_kernel = FALSE
    )

    expect_true(rMVPA:::.fold_cache_enabled())
    expect_true(rMVPA:::.searchlight_geometry_cache_enabled())
    expect_true(rMVPA:::.clustered_nn_fastpath_enabled())
    expect_true(rMVPA:::.matrix_first_roi_enabled())
    expect_true(rMVPA:::.fast_filter_roi_enabled())
    expect_true(rMVPA:::.rsa_fast_kernel_enabled())
    expect_true(rMVPA:::.naive_xdec_fast_kernel_enabled())
  })
})
