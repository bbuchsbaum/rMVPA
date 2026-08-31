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
  expect_true(is.function(env$frperf_ridge_estimators))
  expect_true(is.function(env$frperf_ridge_accuracy_suite))
  expect_true(is.function(env$frperf_ridge_accuracy_summary))
  expect_true(is.function(env$frperf_ridge_accuracy_receipt))
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
      "spectral_block_sse", "center_tcrossprod"
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

test_that("ridge accuracy suite summarizes predeclared linear-fixture margins", {
  script <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "feature_rsa_hotpaths.R"
  )
  skip_if_not(file.exists(script), "benchmark source is not installed")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  rows <- env$frperf_ridge_accuracy_suite(
    seeds = 11:12,
    n_obs = 30L,
    n_features = 5L,
    n_voxels = 8L,
    n_blocks = 5L,
    max_comps = 4L,
    include_glmnet = FALSE
  )
  summary <- env$frperf_ridge_accuracy_summary(rows)
  receipt <- env$frperf_ridge_accuracy_receipt(rows)

  expect_setequal(unique(rows$seed), 11:12)
  expect_true(all(is.finite(rows$holdout_mse)))
  expect_true(all(c(
    "median_mse_ratio", "q90_mse_ratio",
    "median_pattern_delta", "q10_pattern_delta",
    "median_rdm_delta", "q10_rdm_delta",
    "within_linear_fixture_margin"
  ) %in% names(summary)))
  baseline <- summary[summary$implementation == "pls_streaming_loo", ]
  expect_equal(baseline$median_mse_ratio, 1)
  expect_equal(baseline$median_pattern_delta, 0)
  expect_equal(baseline$median_rdm_delta, 0)
  expect_true(baseline$within_linear_fixture_margin)
  expect_equal(nrow(receipt), nrow(summary))
  expect_true(all(c(
    "median_seconds", "median_speedup_vs_pls_loo",
    "source_fingerprint", "worktree_dirty"
  ) %in% names(receipt)))
  expect_true(all(nzchar(receipt$source_fingerprint)))
})

test_that("feature RSA ridge benchmark records speed and held-out accuracy", {
  script <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "feature_rsa_hotpaths.R"
  )
  skip_if_not(file.exists(script), "benchmark source is not installed")
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)

  rows <- env$frperf_ridge_estimators(
    n_obs = 40L,
    n_features = 6L,
    n_voxels = 10L,
    n_blocks = 5L,
    max_comps = 5L,
    iterations = 1L,
    include_glmnet = FALSE
  )

  expect_setequal(
    rows$implementation,
    c(
      "pls_streaming_loo",
      "pls_streaming_blocked",
      "ridge_gcv",
      "ridge_analytic_loo",
      "ridge_blocked"
    )
  )
  expect_true(all(is.finite(rows$elapsed_median_seconds)))
  expect_true(all(rows$elapsed_median_seconds > 0))
  expect_true(all(is.finite(rows$speedup_vs_reference)))
  expect_true(all(is.finite(rows$holdout_mse)))
  expect_true(all(is.finite(rows$holdout_pattern_correlation)))
  expect_true(all(is.finite(rows$holdout_rdm_correlation)))
  ridge_rows <- grepl("^ridge_", rows$implementation)
  expect_true(all(is.finite(rows$selected_lambda[ridge_rows])))
  expect_true(all(is.na(rows$selected_ncomp[ridge_rows])))
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

test_that("feature RSA ridge benchmark receipt is source-bound", {
  source_script <- testthat::test_path(
    "..", "..", "inst", "benchmarks", "feature_rsa_hotpaths.R"
  )
  script <- source_script
  if (!file.exists(script)) {
    script <- system.file(
      "benchmarks", "feature_rsa_hotpaths.R", package = "rMVPA"
    )
  }
  receipt_path <- testthat::test_path(
    "..", "..", "inst", "extdata",
    "feature_rsa_ridge_accuracy_results.csv"
  )
  if (!file.exists(receipt_path)) {
    receipt_path <- system.file(
      "extdata", "feature_rsa_ridge_accuracy_results.csv",
      package = "rMVPA"
    )
  }
  expect_true(file.exists(script))
  expect_true(file.exists(receipt_path))
  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)
  receipt <- utils::read.csv(receipt_path, stringsAsFactors = FALSE)

  expect_setequal(
    receipt$implementation,
    c(
      "pls_streaming_loo", "pls_streaming_blocked", "ridge_gcv",
      "ridge_analytic_loo", "ridge_blocked",
      "glmnet_alpha0_blocked_cv"
    )
  )
  expect_identical(unique(receipt$n_seeds), 20L)
  if (file.exists(source_script)) {
    source_root <- normalizePath(testthat::test_path("..", ".."))
    expect_identical(
      unique(receipt$source_fingerprint),
      env$frperf_source_fingerprint(source_root)
    )
  } else {
    expect_match(
      unique(receipt$source_fingerprint),
      "^[[:xdigit:]]{32}$"
    )
  }
  expect_true(all(receipt[c(
    "omp_threads", "openblas_threads", "mkl_threads", "veclib_threads"
  )] == "1"))
  ridge <- grepl("^ridge_", receipt$implementation)
  expect_true(all(receipt$within_linear_fixture_margin[ridge]))
  expect_true(all(is.finite(receipt$median_seconds)))
  expect_true(all(is.finite(receipt$median_mse_ratio)))
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
  expect_match(
    feature_text,
    "`ncomp_objective = \"mse\"` remains the default",
    fixed = TRUE
  )
  expect_match(feature_text, "pattern_discrimination", fixed = TRUE)
  expect_match(feature_text, "linear-time scorer", fixed = TRUE)
  expect_match(feature_text, "method = \"ridge\"", fixed = TRUE)
  expect_match(feature_text, "Ridge is a\\s+distinct estimator")
  expect_match(feature_text, "feature_rsa_ridge_accuracy_results.csv",
               fixed = TRUE)
  expect_match(feature_text, "profiling did not justify an Rcpp",
               fixed = TRUE)
  expect_match(parallel_text, "Why can Feature RSA still be slow?",
               fixed = TRUE)
  expect_match(parallel_text,
               "does not eliminate private\\s+coefficient")
})
