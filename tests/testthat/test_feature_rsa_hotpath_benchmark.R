test_that("feature RSA benchmark validates parity without timing thresholds", {
  script <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "feature_rsa_hotpaths.R"
  )
  if (!file.exists(script)) {
    script <- system.file(
      "benchmarks", "feature_rsa_hotpaths.R", package = "rMVPA"
    )
  }
  expect_true(file.exists(script))
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  expect_true(is.function(env$frperf_hotpaths))
  expect_true(is.function(env$frperf_worker_payload))
  expect_true(is.function(env$frperf_end_to_end))
  expect_true(is.function(env$frperf_run))

  rows <- env$frperf_hotpaths(
    n_obs = 24L,
    n_features = 5L,
    n_voxels = 8L,
    n_blocks = 4L,
    max_comps = 4L,
    iterations = 1L
  )

  expect_setequal(
    rows$implementation,
    c(
      "compact_matrix_kernel", "streaming_loo", "streaming_blocked",
      "center_tcrossprod"
    )
  )
  expect_true(all(is.finite(rows$elapsed_median_seconds)))
  expect_lte(max(rows$max_abs_prediction_error, na.rm = TRUE), 1e-9)
  expect_lte(max(rows$max_abs_validation_error, na.rm = TRUE), 1e-9)
  selected <- rows[is.finite(rows$selected_ncomp), , drop = FALSE]
  expect_identical(selected$selected_ncomp, selected$reference_selected_ncomp)
  expect_gt(
    rows$payload_reduction_ratio[
      rows$implementation == "compact_matrix_kernel"
    ],
    1
  )
})

test_that("feature RSA worker benchmark isolates design compaction", {
  script <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "feature_rsa_hotpaths.R"
  )
  skip_if_not(file.exists(script), "benchmark source is not installed")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  row <- env$frperf_worker_payload(
    n_obs = 80L,
    n_features = 8L,
    n_blocks = 4L
  )

  expect_identical(row$reference, "dataset_only_removed")
  expect_lt(row$serialized_bytes, row$reference_serialized_bytes)
  expect_gt(row$payload_reduction_ratio, 1.2)
})

test_that("feature RSA documentation distinguishes selector estimands", {
  feature_vignette <- testthat::test_path(
    "..", "..", "vignettes", "Feature_RSA.Rmd"
  )
  parallel_vignette <- testthat::test_path(
    "..", "..", "vignettes", "Parallelism.Rmd"
  )
  skip_if_not(
    file.exists(feature_vignette) && file.exists(parallel_vignette),
    "source vignettes are not installed"
  )
  feature_text <- paste(readLines(feature_vignette, warn = FALSE),
                        collapse = "\n")
  parallel_text <- paste(readLines(parallel_vignette, warn = FALSE),
                         collapse = "\n")

  expect_match(feature_text, "One fit per training block", fixed = TRUE)
  expect_match(feature_text, "One fit per training observation", fixed = TRUE)
  expect_match(feature_text, "equal weight", fixed = TRUE)
  expect_match(feature_text, "fails the fold", fixed = TRUE)
  expect_match(parallel_text, "Why can Feature RSA still be slow?",
               fixed = TRUE)
  expect_match(parallel_text,
               "does not eliminate private\\s+coefficient")
})
