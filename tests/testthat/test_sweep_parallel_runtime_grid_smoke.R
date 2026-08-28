test_that("parallel runtime sweep script writes dry-run manifests", {
  skip_on_cran()

  script_candidate <- testthat::test_path(
    "..", "..", "scripts", "sweep_parallel_runtime_grid.R"
  )
  skip_if_not(
    file.exists(script_candidate),
    "the source-only runtime sweep is not installed with the package"
  )
  script_path <- normalizePath(script_candidate, mustWork = TRUE)
  out_summary <- tempfile("rmvpa-hpc-sweep-summary-", fileext = ".csv")
  out_raw <- tempfile("rmvpa-hpc-sweep-raw-", fileext = ".csv")
  log_dir <- tempfile("rmvpa-hpc-sweep-logs-")

  env <- c(
    "RMVPA_HPC_SWEEP_DRY_RUN=true",
    "RMVPA_HPC_SWEEP_ANALYSES=regional",
    "RMVPA_HPC_SWEEP_BACKENDS=sequential,multisession,callr",
    "RMVPA_HPC_SWEEP_DATA_BACKENDS=default,shard",
    "RMVPA_HPC_SWEEP_LOAD_MODE=source",
    "RMVPA_HPC_SWEEP_WORKER_COUNTS=1,2",
    "RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS=1,2",
    "RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS=1",
    "RMVPA_HPC_SWEEP_BATCH_SIZES=auto,4",
    "RMVPA_HPC_SWEEP_REP=1",
    sprintf("RMVPA_HPC_SWEEP_OUT=%s", out_summary),
    sprintf("RMVPA_HPC_SWEEP_OUT_RAW=%s", out_raw),
    sprintf("RMVPA_HPC_SWEEP_LOG_DIR=%s", log_dir)
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
  expect_equal(status, 0L, info = paste(out, collapse = "\n"))

  expect_true(file.exists(out_summary))
  expect_true(file.exists(out_raw))

  summary_df <- utils::read.csv(out_summary, stringsAsFactors = FALSE)
  raw_df <- utils::read.csv(out_raw, stringsAsFactors = FALSE)

  expect_equal(nrow(raw_df), 48L)
  expect_true(all(c(
    "config_id", "run_id", "analysis", "future_backend", "data_backend",
    "load_mode", "workers", "configured_workers", "omp_threads",
    "blas_threads", "batch_size", "analysis_seconds", "result_signature",
    "frame_bytes_max", "peak_tree_rss_bytes", "status", "message", "log_file"
  ) %in% names(raw_df)))
  expect_true(all(c(
    "config_id", "analysis", "future_backend", "data_backend", "load_mode",
    "workers", "omp_threads", "blas_threads", "batch_size", "n_runs", "n_skip",
    "median_analysis_seconds_ok", "median_frame_bytes_max_ok",
    "median_peak_tree_rss_bytes_ok", "example_log_file"
  ) %in% names(summary_df)))
  expect_true(all(raw_df$analysis == "regional"))
  expect_true(all(raw_df$future_backend %in% c("sequential", "multisession", "callr")))
  expect_true(all(raw_df$data_backend %in% c("default", "shard")))
  expect_true(all(raw_df$load_mode == "source"))
  expect_true(all(raw_df$batch_size %in% c("auto", "4")))
})

test_that("parallel runtime sweep publishes complete child rows atomically", {
  script_candidate <- testthat::test_path(
    "..", "..", "scripts", "sweep_parallel_runtime_grid.R"
  )
  skip_if_not(
    file.exists(script_candidate),
    "the source-only runtime sweep is not installed with the package"
  )
  env <- new.env(parent = globalenv())
  sys.source(script_candidate, envir = env)
  target <- tempfile("rmvpa-sweep-result-", fileext = ".csv")

  expect_false(env$result_file_ready(target))
  expect_true(file.create(target))
  expect_false(env$result_file_ready(target))
  empty <- env$read_result_row(target)
  expect_s3_class(empty, "error")
  expect_match(conditionMessage(empty), "absent or empty")
  unlink(target)

  expected <- data.frame(status = "ok", value = 7L)
  env$atomic_write_csv(expected, target)
  expect_true(env$result_file_ready(target))
  expect_equal(env$read_result_row(target), expected)

  writeLines(c("status,value", "ok,1", "ok,2"), target)
  malformed <- env$read_result_row(target)
  expect_s3_class(malformed, "error")
  expect_match(conditionMessage(malformed), "exactly one")
})
