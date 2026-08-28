# Reproducible rMVPA Future-plan x data-backend characterization.
#
# This driver deliberately runs every configuration in a fresh R process. It
# validates numerical parity before writing the compact receipt used by the
# parallelism vignette. Wall times are descriptive measurements, never package
# test thresholds.
#
# Recommended use from a source checkout after installing an exact tarball into
# an isolated library:
#
#   R_LIBS=/path/to/isolated/library \
#     Rscript inst/benchmarks/parallel_runtime_benchmark.R

parallel_benchmark_required_columns <- function() {
  c(
    "analysis", "future_backend", "data_backend", "load_mode", "workers",
    "configured_workers", "omp_threads", "blas_threads", "batch_size",
    "repetition", "status", "analysis_seconds", "result_keys",
    "result_values", "result_signature", "data_bytes", "frame_bytes_sum",
    "frame_bytes_max", "roi_extract_seconds", "run_future_seconds",
    "n_batches", "package_version", "package_path", "future_version",
    "furrr_version", "shard_version", "r_version", "platform"
  )
}

parallel_benchmark_parse_values <- function(x) {
  if (is.na(x) || !nzchar(x)) return(numeric())
  suppressWarnings(as.numeric(strsplit(x, ";", fixed = TRUE)[[1L]]))
}

parallel_benchmark_hardware <- function() {
  chip <- memory <- NA_character_
  cores <- NA_integer_
  if (identical(Sys.info()[["sysname"]], "Darwin") &&
      nzchar(Sys.which("system_profiler"))) {
    lines <- tryCatch(
      system2("system_profiler", "SPHardwareDataType", stdout = TRUE,
              stderr = FALSE),
      error = function(e) character()
    )
    chip_line <- grep("^[[:space:]]+Chip:", lines, value = TRUE)
    memory_line <- grep("^[[:space:]]+Memory:", lines, value = TRUE)
    cores_line <- grep("^[[:space:]]+Total Number of Cores:", lines,
                       value = TRUE)
    if (length(chip_line)) {
      chip <- sub(".*Chip:[[:space:]]*", "", chip_line[[1L]])
    }
    if (length(memory_line)) {
      memory <- sub(".*Memory:[[:space:]]*", "", memory_line[[1L]])
    }
    if (length(cores_line)) {
      cores <- suppressWarnings(as.integer(sub(
        ".*Total Number of Cores:[[:space:]]*([0-9]+).*", "\\1",
        cores_line[[1L]]
      )))
    }
  }
  list(chip = chip, memory = memory, cores = cores)
}

parallel_benchmark_git_value <- function(args, default = NA_character_) {
  out <- tryCatch(
    system2("git", args, stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  if (length(out)) out[[1L]] else default
}

parallel_benchmark_source_fingerprint <- function(source_root = ".") {
  paths <- file.path(source_root, c(
    "R/mvpa_iterate.R",
    "R/regional.R",
    "R/shard_backend.R",
    "scripts/sweep_parallel_runtime_grid.R",
    "inst/benchmarks/parallel_runtime_benchmark.R"
  ))
  if (!all(file.exists(paths))) return(NA_character_)
  manifest <- tempfile("rmvpa-parallel-source-", fileext = ".txt")
  on.exit(unlink(manifest, force = TRUE), add = TRUE)
  writeLines(
    paste(basename(paths), unname(tools::md5sum(paths)), sep = "="),
    manifest
  )
  unname(as.character(tools::md5sum(manifest)))
}

parallel_benchmark_validate <- function(raw) {
  missing <- setdiff(parallel_benchmark_required_columns(), names(raw))
  if (length(missing)) {
    stop("Raw sweep is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (nrow(raw) < 1L || any(raw$status != "ok")) {
    bad <- unique(raw$status[raw$status != "ok"])
    stop("Every benchmark run must finish with status='ok'; observed: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  if (!all(raw$analysis == "regional") || !all(raw$load_mode == "installed")) {
    stop("The publication receipt requires the installed regional workload.",
         call. = FALSE)
  }
  if (!all(raw$omp_threads == 1L) || !all(raw$blas_threads == 1L)) {
    stop("Process/data-backend comparisons require OMP and BLAS threads = 1.",
         call. = FALSE)
  }
  if (any(raw$n_batches < 2L)) {
    stop("Benchmark did not create multiple batches; Future workers were not exercised.",
         call. = FALSE)
  }
  nonsequential <- raw$future_backend != "sequential"
  if (any(raw$configured_workers[nonsequential] != raw$workers[nonsequential])) {
    stop("A non-sequential plan did not configure the requested worker count.",
         call. = FALSE)
  }
  if (any(raw$configured_workers[!nonsequential] != 1L)) {
    stop("Sequential runs must report one configured worker.", call. = FALSE)
  }

  errors <- numeric()
  for (rep_id in unique(raw$repetition)) {
    rows <- raw[raw$repetition == rep_id, , drop = FALSE]
    if (length(unique(rows$result_keys)) != 1L) {
      stop("Result keys differ within repetition ", rep_id, ".", call. = FALSE)
    }
    reference <- parallel_benchmark_parse_values(rows$result_values[[1L]])
    if (!length(reference) || any(!is.finite(reference))) {
      stop("Reference result is empty or non-finite in repetition ", rep_id,
           ".", call. = FALSE)
    }
    rep_errors <- vapply(rows$result_values, function(value) {
      candidate <- parallel_benchmark_parse_values(value)
      if (length(candidate) != length(reference) || any(!is.finite(candidate))) {
        return(Inf)
      }
      max(abs(candidate - reference))
    }, numeric(1))
    errors <- c(errors, rep_errors)
    if (length(unique(rows$result_signature)) != 1L) {
      stop("Result signatures differ within repetition ", rep_id, ".",
           call. = FALSE)
    }
  }
  max_error <- max(errors)
  if (!is.finite(max_error) || max_error > 1e-12) {
    stop(sprintf("Maximum backend result difference was %.17g.", max_error),
         call. = FALSE)
  }
  max_error
}

parallel_benchmark_summarise <- function(raw_file,
                                         receipt_file,
                                         source_root = ".") {
  raw <- utils::read.csv(raw_file, stringsAsFactors = FALSE,
                         check.names = FALSE)
  max_error <- parallel_benchmark_validate(raw)
  split_key <- interaction(raw$future_backend, raw$data_backend, drop = TRUE)
  groups <- split(raw, split_key)
  order_backends <- c(
    "sequential", "multisession", "multicore", "mirai_multisession"
  )
  hardware <- parallel_benchmark_hardware()
  git_status <- tryCatch(
    system2("git", c("status", "--porcelain"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )

  receipt <- do.call(rbind, lapply(groups, function(rows) {
    data.frame(
      future_backend = rows$future_backend[[1L]],
      data_backend = rows$data_backend[[1L]],
      repetitions = nrow(rows),
      requested_workers = rows$workers[[1L]],
      configured_workers = rows$configured_workers[[1L]],
      omp_threads = rows$omp_threads[[1L]],
      blas_threads = rows$blas_threads[[1L]],
      batch_size = rows$batch_size[[1L]],
      n_batches = rows$n_batches[[1L]],
      analysis_median_seconds = stats::median(rows$analysis_seconds),
      analysis_min_seconds = min(rows$analysis_seconds),
      analysis_max_seconds = max(rows$analysis_seconds),
      frame_bytes_sum_median = stats::median(rows$frame_bytes_sum),
      frame_bytes_max_median = stats::median(rows$frame_bytes_max),
      roi_extract_median_seconds = stats::median(rows$roi_extract_seconds),
      future_map_median_seconds = stats::median(rows$run_future_seconds),
      max_abs_result_error = max_error,
      fixture_data_bytes = stats::median(rows$data_bytes),
      package_version = rows$package_version[[1L]],
      load_mode = rows$load_mode[[1L]],
      future_version = rows$future_version[[1L]],
      furrr_version = rows$furrr_version[[1L]],
      shard_version = rows$shard_version[[1L]],
      r_version = rows$r_version[[1L]],
      platform = rows$platform[[1L]],
      hardware_chip = hardware$chip,
      hardware_memory = hardware$memory,
      hardware_cpu_cores = hardware$cores,
      source_git_head = parallel_benchmark_git_value(c("rev-parse", "HEAD")),
      source_worktree_dirty = length(git_status) > 0L,
      source_fingerprint = parallel_benchmark_source_fingerprint(source_root),
      benchmark_date = as.character(Sys.Date()),
      stringsAsFactors = FALSE
    )
  }))
  receipt <- receipt[
    order(match(receipt$future_backend, order_backends), receipt$data_backend),
    , drop = FALSE
  ]

  default_sum <- setNames(
    receipt$frame_bytes_sum_median[receipt$data_backend == "default"],
    receipt$future_backend[receipt$data_backend == "default"]
  )
  receipt$default_to_shard_frame_sum_ratio <- ifelse(
    receipt$data_backend == "shard",
    default_sum[receipt$future_backend] / receipt$frame_bytes_sum_median,
    1
  )
  receipt$rss_interpretation <- paste(
    "not published: summed process RSS can miss detached workers and",
    "double-count shared or copy-on-write pages"
  )

  dir.create(dirname(receipt_file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(receipt, receipt_file, row.names = FALSE)
  receipt
}

parallel_benchmark_run <- function(
    receipt_file = file.path("inst", "extdata",
                             "parallel_runtime_benchmark_results.csv"),
    repetitions = 3L,
    source_root = ".") {
  sweep_script <- file.path(source_root, "scripts",
                            "sweep_parallel_runtime_grid.R")
  if (!file.exists(sweep_script)) {
    stop("Run this benchmark from the rMVPA source root.", call. = FALSE)
  }
  if (!identical(normalizePath(system.file(package = "rMVPA"), mustWork = TRUE),
                 normalizePath(file.path(.libPaths()[[1L]], "rMVPA"),
                               mustWork = TRUE))) {
    stop(
      "rMVPA must be installed in the first library path; use an isolated R_LIBS.",
      call. = FALSE
    )
  }

  diagnostic_root <- if (dir.exists("/tmp")) "/tmp" else tempdir()
  run_dir <- tempfile(
    "rmvpa-parallel-runtime-",
    tmpdir = diagnostic_root
  )
  dir.create(run_dir, recursive = TRUE)
  completed <- FALSE
  on.exit({
    if (isTRUE(completed)) {
      unlink(run_dir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  raw_file <- file.path(run_dir, "raw.csv")
  summary_file <- file.path(run_dir, "summary.csv")
  log_dir <- file.path(run_dir, "logs")
  env <- c(
    sprintf("R_LIBS=%s", paste(.libPaths(), collapse = .Platform$path.sep)),
    "R_PARALLELLY_MAXWORKERS_LOCALHOST=2",
    "RMVPA_HPC_SWEEP_ANALYSES=regional",
    paste0("RMVPA_HPC_SWEEP_BACKENDS=",
           "sequential,multisession,multicore,mirai_multisession"),
    "RMVPA_HPC_SWEEP_DATA_BACKENDS=default,shard",
    "RMVPA_HPC_SWEEP_LOAD_MODE=installed",
    "RMVPA_HPC_SWEEP_WORKER_COUNTS=2",
    "RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS=1",
    "RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS=1",
    "RMVPA_HPC_SWEEP_BATCH_SIZES=4",
    sprintf("RMVPA_HPC_SWEEP_REP=%d", as.integer(repetitions)),
    "RMVPA_HPC_SWEEP_TIMEOUT_SECONDS=180",
    "RMVPA_HPC_SWEEP_DIMS=28,28,28",
    "RMVPA_HPC_SWEEP_NOBS=72",
    "RMVPA_HPC_SWEEP_NLEVELS=3",
    "RMVPA_HPC_SWEEP_BLOCKS=6",
    "RMVPA_HPC_SWEEP_N_REGIONS=18",
    sprintf("RMVPA_HPC_SWEEP_OUT=%s", summary_file),
    sprintf("RMVPA_HPC_SWEEP_OUT_RAW=%s", raw_file),
    sprintf("RMVPA_HPC_SWEEP_LOG_DIR=%s", log_dir)
  )
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    sweep_script,
    env = env
  )
  if (!identical(as.integer(status), 0L)) {
    stop(
      "Parallel runtime sweep failed with status ", status,
      "; diagnostic files retained at ", run_dir, ".",
      call. = FALSE
    )
  }
  result <- parallel_benchmark_summarise(
    raw_file, receipt_file, source_root
  )
  completed <- TRUE
  result
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  receipt_file <- if (length(args)) args[[1L]] else file.path(
    "inst", "extdata", "parallel_runtime_benchmark_results.csv"
  )
  result <- parallel_benchmark_run(receipt_file = receipt_file)
  print(result[, c(
    "future_backend", "data_backend", "analysis_median_seconds",
    "frame_bytes_sum_median", "max_abs_result_error"
  )], row.names = FALSE)
}
