test_that("parallel runtime sweep script writes dry-run manifests", {
  skip_on_cran()

  script_path <- normalizePath(
    testthat::test_path("..", "..", "scripts", "sweep_parallel_runtime_grid.R"),
    mustWork = TRUE
  )
  out_summary <- tempfile("rmvpa-hpc-sweep-summary-", fileext = ".csv")
  out_raw <- tempfile("rmvpa-hpc-sweep-raw-", fileext = ".csv")
  log_dir <- tempfile("rmvpa-hpc-sweep-logs-")

  env <- c(
    "RMVPA_HPC_SWEEP_DRY_RUN=true",
    "RMVPA_HPC_SWEEP_ANALYSES=regional",
    "RMVPA_HPC_SWEEP_BACKENDS=sequential,multisession,callr",
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

  expect_equal(nrow(raw_df), 24L)
  expect_true(all(c(
    "config_id", "run_id", "analysis", "future_backend", "workers", "omp_threads",
    "blas_threads", "batch_size", "status", "message", "log_file"
  ) %in% names(raw_df)))
  expect_true(all(c(
    "config_id", "analysis", "future_backend", "workers", "omp_threads",
    "blas_threads", "batch_size", "n_runs", "n_skip", "example_log_file"
  ) %in% names(summary_df)))
  expect_true(all(raw_df$analysis == "regional"))
  expect_true(all(raw_df$future_backend %in% c("sequential", "multisession", "callr")))
  expect_true(all(raw_df$batch_size %in% c("auto", "4")))
})
