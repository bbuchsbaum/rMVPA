test_that("feature RSA hot paths retain their measured speed advantage", {
  skip_on_cran()
  skip_if_not_perf_tests()

  script <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "feature_rsa_hotpaths.R"
  )
  skip_if_not(file.exists(script), "benchmark source is not installed")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  rows <- env$frperf_hotpaths(iterations = 3L)
  guarded <- rows$implementation %in% c(
    "compact_matrix_kernel", "streaming_loo", "center_tcrossprod"
  )
  message(paste(
    sprintf(
      "%s %.2fx",
      rows$implementation[guarded],
      rows$speedup_vs_reference[guarded]
    ),
    collapse = "; "
  ))

  # A broad regression guard, enabled only in the dedicated performance lane.
  # Scientific parity is enforced independently in ordinary tests.
  expect_true(all(
    rows$elapsed_median_seconds[guarded] <=
      1.25 * rows$reference_median_seconds[guarded]
  ))
})
