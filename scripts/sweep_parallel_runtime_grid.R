suppressPackageStartupMessages({
  library(parallel)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

flag_enabled <- function(name, default = FALSE) {
  raw <- trimws(Sys.getenv(name, if (isTRUE(default)) "true" else "false"))
  tolower(raw) %in% c("1", "true", "yes", "on")
}

env_int <- function(name, default) {
  val <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(val) || val < 1L) default else val
}

env_csv <- function(name, default = character()) {
  raw <- trimws(Sys.getenv(name, ""))
  if (!nzchar(raw)) {
    return(default)
  }
  vals <- trimws(unlist(strsplit(raw, ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

env_choice <- function(name, choices, default) {
  val <- tolower(trimws(Sys.getenv(name, default)))
  if (!nzchar(val)) {
    val <- default
  }
  if (!val %in% choices) {
    stop(sprintf(
      "Invalid value for %s: %s (expected one of %s)",
      name, val, paste(choices, collapse = ", ")
    ), call. = FALSE)
  }
  val
}

parse_positive_ints <- function(values, what) {
  out <- suppressWarnings(as.integer(values))
  if (length(out) == 0L || any(is.na(out)) || any(out < 1L)) {
    stop(sprintf("%s must contain only positive integers.", what), call. = FALSE)
  }
  unique(out)
}

parse_batch_values <- function(values) {
  vals <- trimws(values)
  if (length(vals) == 0L) {
    return("auto")
  }
  valid <- vapply(vals, function(x) {
    identical(x, "auto") || (!is.na(suppressWarnings(as.integer(x))) &&
      suppressWarnings(as.integer(x)) >= 1L)
  }, logical(1))
  if (!all(valid)) {
    stop("RMVPA_HPC_SWEEP_BATCH_SIZE must contain 'auto' and/or positive integers.", call. = FALSE)
  }
  unique(vals)
}

current_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) < 1L) {
    stop("Unable to determine script path from commandArgs().", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
}

safe_slug <- function(x) {
  out <- gsub("[^A-Za-z0-9._-]+", "-", as.character(x))
  out <- gsub("-+", "-", out)
  out <- sub("^-", "", out)
  out <- sub("-$", "", out)
  if (!nzchar(out)) "case" else out
}

timestamp_iso <- function(x = Sys.time()) {
  format(as.POSIXct(x, tz = ""), "%Y-%m-%dT%H:%M:%S%z")
}

is_process_running <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (!is.finite(pid) || is.na(pid) || pid <= 0L) {
    return(FALSE)
  }
  status <- suppressWarnings(
    system2("kill", c("-0", as.character(pid)), stdout = FALSE, stderr = FALSE)
  )
  identical(as.integer(status), 0L)
}

terminate_process <- function(pid, grace_seconds = 2) {
  pid <- suppressWarnings(as.integer(pid))
  if (!is.finite(pid) || is.na(pid) || pid <= 0L) {
    return(invisible(FALSE))
  }
  if (!is_process_running(pid)) {
    return(invisible(TRUE))
  }
  suppressWarnings(tools::pskill(pid, tools::SIGTERM))
  deadline <- Sys.time() + grace_seconds
  while (is_process_running(pid) && Sys.time() < deadline) {
    Sys.sleep(0.1)
  }
  if (is_process_running(pid)) {
    suppressWarnings(tools::pskill(pid, tools::SIGKILL))
  }
  invisible(!is_process_running(pid))
}

tail_log <- function(path, n = 40L) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0L) {
    return("")
  }
  paste(utils::tail(lines, n), collapse = "\n")
}

backend_availability <- function(backend) {
  backend <- tolower(backend)
  if (identical(backend, "sequential")) {
    return(list(ok = TRUE, reason = NA_character_))
  }
  if (identical(backend, "multisession")) {
    return(list(ok = TRUE, reason = NA_character_))
  }
  if (identical(backend, "multicore")) {
    if (!isTRUE(future::supportsMulticore())) {
      return(list(ok = FALSE, reason = "future::supportsMulticore() is FALSE"))
    }
    return(list(ok = TRUE, reason = NA_character_))
  }
  if (identical(backend, "callr")) {
    if (!requireNamespace("future.callr", quietly = TRUE)) {
      return(list(ok = FALSE, reason = "future.callr is not installed"))
    }
    return(list(ok = TRUE, reason = NA_character_))
  }
  if (identical(backend, "mirai_multisession")) {
    if (!requireNamespace("future.mirai", quietly = TRUE)) {
      return(list(ok = FALSE, reason = "future.mirai is not installed"))
    }
    if (!exists("mirai_multisession", envir = asNamespace("future.mirai"), inherits = FALSE)) {
      return(list(ok = FALSE, reason = "future.mirai::mirai_multisession() is unavailable"))
    }
    return(list(ok = TRUE, reason = NA_character_))
  }
  list(ok = FALSE, reason = sprintf("Unknown backend '%s'", backend))
}

configure_future_plan <- function(backend, workers) {
  backend <- tolower(backend)
  if (identical(backend, "sequential")) {
    future::plan(future::sequential)
    return("sequential")
  }
  if (identical(backend, "multisession")) {
    future::plan(future::multisession, workers = workers)
    return("multisession")
  }
  if (identical(backend, "multicore")) {
    if (!isTRUE(future::supportsMulticore())) {
      stop("future::multicore is not supported in this environment.", call. = FALSE)
    }
    future::plan(future::multicore, workers = workers)
    return("multicore")
  }
  if (identical(backend, "callr")) {
    if (!requireNamespace("future.callr", quietly = TRUE)) {
      stop("future.callr is not installed.", call. = FALSE)
    }
    future::plan(get("callr", envir = asNamespace("future.callr")), workers = workers)
    return("callr")
  }
  if (identical(backend, "mirai_multisession")) {
    if (!requireNamespace("future.mirai", quietly = TRUE)) {
      stop("future.mirai is not installed.", call. = FALSE)
    }
    plan_fun <- get("mirai_multisession", envir = asNamespace("future.mirai"))
    future::plan(plan_fun, workers = workers)
    return("mirai_multisession")
  }
  stop(sprintf("Unsupported backend '%s'", backend), call. = FALSE)
}

build_mvpa_spec <- function() {
  ds <- gen_sample_dataset(D = c(7, 7, 7), nobs = 48, nlevels = 3, blocks = 4)
  cv <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
  list(dataset = ds, model_spec = mspec)
}

build_region_mask <- function(dataset, n_regions = 18L) {
  mask_template <- dataset$dataset$mask
  mask_vec <- integer(length(mask_template))
  inside <- which(as.vector(mask_template) > 0)
  mask_vec[inside] <- sample.int(n_regions, length(inside), replace = TRUE)
  neuroim2::NeuroVol(mask_vec, neuroim2::space(mask_template))
}

run_probe_case <- function(config) {
  suppressPackageStartupMessages({
    library(pkgload)
  })

  pkgload::load_all(".", quiet = TRUE)

  if (requireNamespace("futile.logger", quietly = TRUE)) {
    futile.logger::flog.threshold(futile.logger::ERROR)
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  configured_backend <- configure_future_plan(config$future_backend, config$workers)

  fixture <- build_mvpa_spec()
  batch_size <- if (identical(config$batch_size, "auto")) NULL else as.integer(config$batch_size)

  if (identical(config$analysis, "regional")) {
    region_mask <- build_region_mask(fixture$dataset, n_regions = 18L)
    args <- list(
      model_spec = fixture$model_spec,
      region_mask = region_mask,
      verbose = FALSE,
      backend = "default"
    )
    if (!is.null(batch_size)) {
      args$batch_size <- batch_size
    }
    res <- do.call(run_regional, args)
    perf_rows <- nrow(res$performance_table %||% tibble::tibble())
    prediction_rows <- nrow(res$prediction_table %||% tibble::tibble())
    list(
      result_class = paste(class(res), collapse = ";"),
      perf_rows = perf_rows,
      prediction_rows = prediction_rows,
      fit_count = length(res$fits %||% list()),
      configured_backend = configured_backend
    )
  } else if (identical(config$analysis, "searchlight")) {
    args <- list(
      model_spec = fixture$model_spec,
      radius = config$radius,
      method = "standard",
      verbose = FALSE,
      backend = "default"
    )
    if (!is.null(batch_size)) {
      args$batch_size <- batch_size
    }
    res <- do.call(run_searchlight, args)
    list(
      result_class = paste(class(res), collapse = ";"),
      perf_rows = length(res),
      prediction_rows = NA_integer_,
      fit_count = NA_integer_,
      configured_backend = configured_backend
    )
  } else {
    stop(sprintf("Unsupported analysis '%s'", config$analysis), call. = FALSE)
  }
}

config_from_env <- function() {
  list(
    analysis = env_choice("RMVPA_HPC_SWEEP_ANALYSIS", c("regional", "searchlight"), "regional"),
    future_backend = env_choice(
      "RMVPA_HPC_SWEEP_FUTURE_BACKEND",
      c("sequential", "multisession", "multicore", "callr", "mirai_multisession"),
      "sequential"
    ),
    workers = env_int("RMVPA_HPC_SWEEP_WORKERS", 1L),
    omp_threads = env_int("RMVPA_HPC_SWEEP_OMP_THREADS", 1L),
    blas_threads = env_int("RMVPA_HPC_SWEEP_BLAS_THREADS", 1L),
    batch_size = Sys.getenv("RMVPA_HPC_SWEEP_BATCH_SIZE", "auto"),
    repetition = env_int("RMVPA_HPC_SWEEP_REPETITION", 1L),
    timeout_seconds = env_int("RMVPA_HPC_SWEEP_TIMEOUT_SECONDS", 120L),
    radius = env_int("RMVPA_HPC_SWEEP_RADIUS", 3L),
    log_file = Sys.getenv("RMVPA_HPC_SWEEP_LOG_FILE", ""),
    result_file = Sys.getenv("RMVPA_HPC_SWEEP_RESULT_FILE", "")
  )
}

config_id <- function(config) {
  safe_slug(sprintf(
    "%s-%s-w%d-omp%d-blas%d-batch%s",
    config$analysis,
    config$future_backend,
    config$workers,
    config$omp_threads,
    config$blas_threads,
    config$batch_size
  ))
}

run_id <- function(config) {
  safe_slug(sprintf("%s-r%d", config_id(config), config$repetition))
}

child_main <- function() {
  config <- config_from_env()
  started_at <- Sys.time()
  pid <- Sys.getpid()
  host <- Sys.info()[["nodename"]] %||% NA_character_
  status <- "ok"
  message_text <- ""
  perf_rows <- NA_integer_
  prediction_rows <- NA_integer_
  fit_count <- NA_integer_
  result_class <- NA_character_
  configured_backend <- NA_character_

  result <- tryCatch({
    set.seed(9000L + config$repetition)
    run_probe_case(config)
  }, error = function(e) {
    status <<- "error"
    message_text <<- conditionMessage(e)
    NULL
  })

  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  if (!is.null(result)) {
    perf_rows <- as.integer(result$perf_rows %||% NA_integer_)
    prediction_rows <- as.integer(result$prediction_rows %||% NA_integer_)
    fit_count <- as.integer(result$fit_count %||% NA_integer_)
    result_class <- as.character(result$result_class %||% NA_character_)
    configured_backend <- as.character(result$configured_backend %||% NA_character_)
  }

  row <- data.frame(
    config_id = config_id(config),
    run_id = run_id(config),
    analysis = config$analysis,
    future_backend = config$future_backend,
    configured_backend = configured_backend,
    workers = as.integer(config$workers),
    omp_threads = as.integer(config$omp_threads),
    blas_threads = as.integer(config$blas_threads),
    batch_size = as.character(config$batch_size),
    repetition = as.integer(config$repetition),
    timeout_seconds = as.integer(config$timeout_seconds),
    status = status,
    elapsed_seconds = elapsed,
    perf_rows = perf_rows,
    prediction_rows = prediction_rows,
    fit_count = fit_count,
    result_class = result_class,
    message = message_text,
    pid = as.integer(pid),
    host = host,
    started_at = timestamp_iso(started_at),
    finished_at = timestamp_iso(Sys.time()),
    log_file = config$log_file,
    stringsAsFactors = FALSE
  )

  if (nzchar(config$result_file)) {
    utils::write.csv(row, config$result_file, row.names = FALSE)
  } else {
    print(row)
  }
  invisible(row)
}

expand_configs <- function() {
  analyses <- env_csv("RMVPA_HPC_SWEEP_ANALYSES", default = "regional")
  backends <- env_csv("RMVPA_HPC_SWEEP_BACKENDS", default = c("sequential", "multisession", "multicore"))
  workers <- parse_positive_ints(env_csv("RMVPA_HPC_SWEEP_WORKER_COUNTS", default = c("1", "2")), "RMVPA_HPC_SWEEP_WORKER_COUNTS")
  omp_threads <- parse_positive_ints(env_csv("RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS", default = c("1", "2")), "RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS")
  blas_threads <- parse_positive_ints(env_csv("RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS", default = c("1")), "RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS")
  batch_sizes <- parse_batch_values(env_csv("RMVPA_HPC_SWEEP_BATCH_SIZES", default = "auto"))
  repetitions <- seq_len(env_int("RMVPA_HPC_SWEEP_REP", 1L))
  timeout_seconds <- env_int("RMVPA_HPC_SWEEP_TIMEOUT_SECONDS", 120L)
  radius <- env_int("RMVPA_HPC_SWEEP_RADIUS", 3L)

  grid <- expand.grid(
    analysis = analyses,
    future_backend = backends,
    workers = workers,
    omp_threads = omp_threads,
    blas_threads = blas_threads,
    batch_size = batch_sizes,
    repetition = repetitions,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )
  grid$timeout_seconds <- timeout_seconds
  grid$radius <- radius
  grid
}

launch_child <- function(script_path, config, log_file, result_file) {
  env <- c(
    "RMVPA_HPC_SWEEP_MODE=child",
    sprintf("RMVPA_HPC_SWEEP_ANALYSIS=%s", config$analysis),
    sprintf("RMVPA_HPC_SWEEP_FUTURE_BACKEND=%s", config$future_backend),
    sprintf("RMVPA_HPC_SWEEP_WORKERS=%d", as.integer(config$workers)),
    sprintf("RMVPA_HPC_SWEEP_OMP_THREADS=%d", as.integer(config$omp_threads)),
    sprintf("RMVPA_HPC_SWEEP_BLAS_THREADS=%d", as.integer(config$blas_threads)),
    sprintf("RMVPA_HPC_SWEEP_BATCH_SIZE=%s", config$batch_size),
    sprintf("RMVPA_HPC_SWEEP_REPETITION=%d", as.integer(config$repetition)),
    sprintf("RMVPA_HPC_SWEEP_TIMEOUT_SECONDS=%d", as.integer(config$timeout_seconds)),
    sprintf("RMVPA_HPC_SWEEP_RADIUS=%d", as.integer(config$radius)),
    sprintf("RMVPA_HPC_SWEEP_LOG_FILE=%s", log_file),
    sprintf("RMVPA_HPC_SWEEP_RESULT_FILE=%s", result_file),
    sprintf("OMP_NUM_THREADS=%d", as.integer(config$omp_threads)),
    sprintf("OPENBLAS_NUM_THREADS=%d", as.integer(config$blas_threads)),
    sprintf("MKL_NUM_THREADS=%d", as.integer(config$blas_threads)),
    sprintf("BLIS_NUM_THREADS=%d", as.integer(config$blas_threads)),
    sprintf("VECLIB_MAXIMUM_THREADS=%d", as.integer(config$blas_threads)),
    "R_PROGRESSR_ENABLE=false"
  )
  pid_file <- tempfile("rmvpa-hpc-pid-", tmpdir = dirname(log_file), fileext = ".txt")
  launcher_file <- tempfile("rmvpa-hpc-launch-", tmpdir = dirname(log_file), fileext = ".sh")
  writeLines(
    c(
      "#!/bin/sh",
      sprintf(
        "%s --vanilla %s > %s 2>&1 &",
        shQuote(file.path(R.home("bin"), "Rscript")),
        shQuote(script_path),
        shQuote(log_file)
      ),
      sprintf("echo $! > %s", shQuote(pid_file))
    ),
    con = launcher_file
  )
  on.exit(unlink(c(pid_file, launcher_file), force = TRUE), add = TRUE)
  status <- suppressWarnings(system2("sh", launcher_file, stdout = FALSE, stderr = FALSE, env = env))
  if (!identical(as.integer(status), 0L) || !file.exists(pid_file)) {
    stop(sprintf("Failed to launch child wrapper for %s", config_id(config)), call. = FALSE)
  }
  pid_text <- readLines(pid_file, warn = FALSE)
  pid <- suppressWarnings(as.integer(pid_text[[1L]] %||% NA_character_))
  if (!is.finite(pid) || is.na(pid) || pid <= 0L) {
    stop(sprintf("Failed to launch child process for %s", config_id(config)), call. = FALSE)
  }
  pid
}

wait_for_result <- function(pid, result_file, timeout_seconds) {
  started <- Sys.time()
  repeat {
    if (file.exists(result_file)) {
      return("finished")
    }
    if (!is_process_running(pid)) {
      return("exited")
    }
    if (as.numeric(difftime(Sys.time(), started, units = "secs")) > timeout_seconds) {
      return("timeout")
    }
    Sys.sleep(0.25)
  }
}

row_from_config <- function(config, status, log_file, message = NA_character_) {
  data.frame(
    config_id = config_id(config),
    run_id = run_id(config),
    analysis = config$analysis,
    future_backend = config$future_backend,
    configured_backend = NA_character_,
    workers = as.integer(config$workers),
    omp_threads = as.integer(config$omp_threads),
    blas_threads = as.integer(config$blas_threads),
    batch_size = as.character(config$batch_size),
    repetition = as.integer(config$repetition),
    timeout_seconds = as.integer(config$timeout_seconds),
    status = status,
    elapsed_seconds = NA_real_,
    perf_rows = NA_integer_,
    prediction_rows = NA_integer_,
    fit_count = NA_integer_,
    result_class = NA_character_,
    message = message,
    pid = NA_integer_,
    host = NA_character_,
    started_at = NA_character_,
    finished_at = NA_character_,
    log_file = log_file,
    stringsAsFactors = FALSE
  )
}

summarise_rows <- function(rows) {
  split_rows <- split(rows, rows$config_id)
  out <- lapply(split_rows, function(df) {
    ok <- df$status == "ok"
    problem_idx <- which(!is.na(df$message) & nzchar(df$message))
    first_problem <- if (length(problem_idx) > 0L) df$message[[problem_idx[[1L]]]] else NA_character_
    data.frame(
      config_id = df$config_id[[1L]],
      analysis = df$analysis[[1L]],
      future_backend = df$future_backend[[1L]],
      workers = df$workers[[1L]],
      omp_threads = df$omp_threads[[1L]],
      blas_threads = df$blas_threads[[1L]],
      batch_size = df$batch_size[[1L]],
      n_runs = nrow(df),
      n_ok = sum(df$status == "ok"),
      n_error = sum(df$status == "error"),
      n_timeout = sum(df$status == "timeout"),
      n_crash = sum(df$status == "crash"),
      n_skip = sum(df$status == "skip"),
      median_elapsed_ok = if (any(ok)) stats::median(df$elapsed_seconds[ok], na.rm = TRUE) else NA_real_,
      example_message = first_problem,
      example_log_file = df$log_file[[1L]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

main <- function() {
  mode <- env_choice("RMVPA_HPC_SWEEP_MODE", c("parent", "child"), "parent")
  if (identical(mode, "child")) {
    child_main()
    return(invisible(NULL))
  }

  if (.Platform$OS.type == "windows") {
    stop("This sweep script currently supports Unix-like systems only.", call. = FALSE)
  }

  script_path <- current_script_path()
  raw_out <- Sys.getenv("RMVPA_HPC_SWEEP_OUT_RAW", file.path(getwd(), ".tmp", "hpc_parallel_sweep_raw.csv"))
  summary_out <- Sys.getenv("RMVPA_HPC_SWEEP_OUT", file.path(getwd(), ".tmp", "hpc_parallel_sweep_summary.csv"))
  log_dir <- Sys.getenv("RMVPA_HPC_SWEEP_LOG_DIR", file.path(getwd(), ".tmp", "hpc_parallel_sweep_logs"))
  dry_run <- flag_enabled("RMVPA_HPC_SWEEP_DRY_RUN", default = FALSE)

  dir.create(dirname(raw_out), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(summary_out), recursive = TRUE, showWarnings = FALSE)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

  configs <- expand_configs()
  if (nrow(configs) < 1L) {
    stop("No configurations selected for sweep.", call. = FALSE)
  }

  message(sprintf("Parallel sweep manifest has %d configuration(s).", nrow(configs)))

  if (isTRUE(dry_run)) {
    rows <- vector("list", nrow(configs))
    for (i in seq_len(nrow(configs))) {
      cfg <- as.list(configs[i, , drop = FALSE])
      cfg$batch_size <- as.character(cfg$batch_size)
      log_file <- file.path(log_dir, sprintf("%03d-%s.log", i, run_id(cfg)))
      avail <- backend_availability(cfg$future_backend)
      rows[[i]] <- row_from_config(
        cfg,
        status = if (isTRUE(avail$ok)) "planned" else "skip",
        log_file = log_file,
        message = avail$reason
      )
    }
    raw_df <- do.call(rbind, rows)
    summary_df <- summarise_rows(transform(raw_df, status = ifelse(status == "planned", "skip", status)))
    utils::write.csv(raw_df, raw_out, row.names = FALSE)
    utils::write.csv(summary_df, summary_out, row.names = FALSE)
    message(sprintf("Dry run manifest written to %s and %s", summary_out, raw_out))
    return(invisible(list(raw = raw_df, summary = summary_df)))
  }

  rows <- vector("list", nrow(configs))
  for (i in seq_len(nrow(configs))) {
    cfg <- as.list(configs[i, , drop = FALSE])
    cfg$batch_size <- as.character(cfg$batch_size)
    log_file <- file.path(log_dir, sprintf("%03d-%s.log", i, run_id(cfg)))
    result_file <- file.path(log_dir, sprintf("%03d-%s.csv", i, run_id(cfg)))
    avail <- backend_availability(cfg$future_backend)

    if (!isTRUE(avail$ok)) {
      rows[[i]] <- row_from_config(cfg, status = "skip", log_file = log_file, message = avail$reason)
      next
    }

    unlink(result_file, force = TRUE)
    message(sprintf(
      "[%d/%d] %s backend=%s workers=%d omp=%d blas=%d batch=%s",
      i, nrow(configs), cfg$analysis, cfg$future_backend, cfg$workers,
      cfg$omp_threads, cfg$blas_threads, cfg$batch_size
    ))

    started <- Sys.time()
    pid <- tryCatch(
      launch_child(script_path, cfg, log_file = log_file, result_file = result_file),
      error = function(e) e
    )
    if (inherits(pid, "error")) {
      rows[[i]] <- row_from_config(cfg, status = "crash", log_file = log_file, message = conditionMessage(pid))
      next
    }

    state <- wait_for_result(pid, result_file, timeout_seconds = cfg$timeout_seconds)

    if (identical(state, "timeout")) {
      terminate_process(pid)
      msg <- sprintf(
        "Timed out after %ds. Log tail:\n%s",
        cfg$timeout_seconds,
        tail_log(log_file, n = 20L)
      )
      row <- row_from_config(cfg, status = "timeout", log_file = log_file, message = msg)
      row$pid <- as.integer(pid)
      row$started_at <- timestamp_iso(started)
      row$finished_at <- timestamp_iso(Sys.time())
      rows[[i]] <- row
      next
    }

    if (file.exists(result_file)) {
      row <- utils::read.csv(result_file, stringsAsFactors = FALSE)
      row$log_file <- log_file
      rows[[i]] <- row
      next
    }

    msg <- sprintf(
      "Child exited without writing a result row. Log tail:\n%s",
      tail_log(log_file, n = 20L)
    )
    row <- row_from_config(cfg, status = "crash", log_file = log_file, message = msg)
    row$pid <- as.integer(pid)
    row$started_at <- timestamp_iso(started)
    row$finished_at <- timestamp_iso(Sys.time())
    rows[[i]] <- row
  }

  raw_df <- do.call(rbind, rows)
  summary_df <- summarise_rows(raw_df)

  utils::write.csv(raw_df, raw_out, row.names = FALSE)
  utils::write.csv(summary_df, summary_out, row.names = FALSE)

  message(sprintf("Wrote summary to %s", summary_out))
  message(sprintf("Wrote raw results to %s", raw_out))

  invisible(list(raw = raw_df, summary = summary_df))
}

main()
