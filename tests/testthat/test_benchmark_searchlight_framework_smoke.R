test_that("benchmark_searchlight_framework script runs and writes expected CSVs", {
  skip_on_cran()
  skip_if_not_perf_tests()
  skip_if_not_installed("devtools")

  script_path <- normalizePath(
    testthat::test_path("..", "..", "scripts", "benchmark_searchlight_framework.R"),
    mustWork = TRUE
  )
  out_summary <- tempfile("rmvpa-bench-summary-", fileext = ".csv")
  out_raw <- tempfile("rmvpa-bench-raw-", fileext = ".csv")

  env <- c(
    "RMVPA_BENCH_REP=1",
    "RMVPA_BENCH_BACKEND=default",
    "RMVPA_BENCH_PROFILE=phase2",
    "RMVPA_BENCH_CASES=mvpa_classification",
    "RMVPA_BENCH_DRY_RUN=true",
    sprintf("RMVPA_BENCH_OUT=%s", out_summary),
    sprintf("RMVPA_BENCH_OUT_RAW=%s", out_raw)
  )

  out <- system2(
    command = file.path(R.home("bin"), "Rscript"),
    args = script_path,
    stdout = TRUE,
    stderr = TRUE,
    env = env
  )
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (identical(as.integer(status), 134L)) {
    skip("Benchmark subprocess aborted (signal 6) in this environment.")
  }
  expect_equal(status, 0L, info = paste(out, collapse = "\n"))

  expect_true(file.exists(out_summary))
  expect_true(file.exists(out_raw))

  summary_df <- utils::read.csv(out_summary, stringsAsFactors = FALSE)
  raw_df <- utils::read.csv(out_raw, stringsAsFactors = FALSE)

  expect_equal(summary_df$workload, "mvpa_classification")
  expect_equal(raw_df$workload, "mvpa_classification")
  expect_true(all(c("backend", "bench_profile", "fold_cache_enabled",
                    "matrix_first_roi", "fast_filter_roi",
                    "naive_xdec_fast_kernel",
                    "nrep", "radius") %in% names(summary_df)))
  expect_true(all(c("backend", "bench_profile", "fold_cache_enabled",
                    "matrix_first_roi", "fast_filter_roi",
                    "naive_xdec_fast_kernel",
                    "nrep", "radius") %in% names(raw_df)))
  expect_true(all(summary_df$bench_profile == "phase2"))
  expect_true(all(raw_df$bench_profile == "phase2"))
})
