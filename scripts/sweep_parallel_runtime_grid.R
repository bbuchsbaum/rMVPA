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

parse_dims <- function(value, what = "RMVPA_HPC_SWEEP_DIMS") {
  raw <- trimws(unlist(strsplit(value, ",", fixed = TRUE)))
  dims <- suppressWarnings(as.integer(raw))
  if (length(dims) != 3L || any(is.na(dims)) || any(dims < 1L)) {
    stop(sprintf("%s must contain exactly three dimensions.", what), call. = FALSE)
  }
  dims
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

result_file_ready <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  info <- file.info(path)
  is.finite(info$size) && info$size > 0L
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(
    paste0(".", basename(path), "-"),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) {
    stop("Failed to publish result file atomically: ", path, call. = FALSE)
  }
  invisible(path)
}

read_result_row <- function(path) {
  tryCatch({
    if (!result_file_ready(path)) {
      stop("result file is absent or empty", call. = FALSE)
    }
    row <- utils::read.csv(path, stringsAsFactors = FALSE)
    if (nrow(row) != 1L) {
      stop("result file must contain exactly one data row", call. = FALSE)
    }
    row
  }, error = function(e) e)
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

data_backend_availability <- function(backend) {
  backend <- tolower(backend)
  if (identical(backend, "default")) {
    return(list(ok = TRUE, reason = NA_character_))
  }
  if (identical(backend, "shard")) {
    if (!requireNamespace("shard", quietly = TRUE)) {
      return(list(ok = FALSE, reason = "shard is not installed"))
    }
    if (.Platform$OS.type == "windows") {
      return(list(ok = FALSE, reason = "shard requires POSIX shared memory"))
    }
    return(list(ok = TRUE, reason = NA_character_))
  }
  list(ok = FALSE, reason = sprintf("Unknown rMVPA data backend '%s'", backend))
}

load_rmvpa <- function(mode) {
  mode <- match.arg(mode, c("source", "installed"))
  if (identical(mode, "source")) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("load_mode='source' requires pkgload.", call. = FALSE)
    }
    pkgload::load_all(".", quiet = TRUE)
  } else {
    suppressPackageStartupMessages(library(rMVPA))
  }
  invisible(mode)
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

build_mvpa_spec <- function(config) {
  ds <- gen_sample_dataset(
    D = config$dims,
    nobs = config$nobs,
    nlevels = config$nlevels,
    blocks = config$blocks
  )
  cv <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(
    model = load_model(config$model),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv,
    return_predictions = FALSE
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

md5_object <- function(x) {
  path <- tempfile("rmvpa-result-signature-", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  saveRDS(x, path, version = 2)
  unname(as.character(tools::md5sum(path)))
}

canonical_result <- function(result, analysis) {
  if (identical(analysis, "regional")) {
    tbl <- as.data.frame(result$performance_table %||% data.frame())
    if (nrow(tbl) < 1L) {
      return(list(keys = character(), values = numeric()))
    }
    id_col <- if ("roinum" %in% names(tbl)) "roinum" else names(tbl)[[1L]]
    tbl <- tbl[order(tbl[[id_col]]), , drop = FALSE]
    metric_cols <- setdiff(names(tbl), id_col)
    metric_cols <- metric_cols[vapply(tbl[metric_cols], is.numeric, logical(1))]
    if (length(metric_cols) < 1L) {
      return(list(keys = character(), values = numeric()))
    }
    mat <- as.matrix(tbl[, metric_cols, drop = FALSE])
    values <- as.numeric(t(mat))
    keys <- unlist(lapply(seq_len(nrow(tbl)), function(i) {
      paste(tbl[[id_col]][[i]], metric_cols, sep = ":")
    }), use.names = FALSE)
    return(list(keys = keys, values = values))
  }

  maps <- result$results %||% list()
  if (length(maps) < 1L) {
    return(list(keys = character(), values = numeric()))
  }
  map_names <- sort(names(maps))
  if (is.null(map_names) || any(!nzchar(map_names))) {
    map_names <- as.character(seq_along(maps))
  } else {
    maps <- maps[map_names]
  }
  keys <- character()
  values <- numeric()
  for (i in seq_along(maps)) {
    vals <- tryCatch(
      as.numeric(neuroim2::values(maps[[i]])),
      error = function(e) suppressWarnings(as.numeric(maps[[i]]))
    )
    keep <- which(is.finite(vals))
    keys <- c(keys, paste(map_names[[i]], keep, sep = ":"))
    values <- c(values, vals[keep])
  }
  list(keys = keys, values = values)
}

profile_receipt <- function(result) {
  timing <- attr(result, "timing", exact = TRUE)
  iterate <- timing$iterate %||% timing
  totals <- iterate$totals %||% list()
  list(
    frame_bytes_sum = as.numeric(totals$frame_bytes_sum %||% NA_real_),
    frame_bytes_max = as.numeric(totals$frame_bytes_max %||% NA_real_),
    roi_extract_seconds = as.numeric(totals$roi_extract_seconds %||% NA_real_),
    run_future_seconds = as.numeric(totals$run_future_seconds %||% NA_real_),
    n_batches = as.integer(iterate$n_batches %||% NA_integer_)
  )
}

run_probe_case <- function(config) {
  load_rmvpa(config$load_mode)

  if (requireNamespace("futile.logger", quietly = TRUE)) {
    futile.logger::flog.threshold(futile.logger::ERROR)
  }

  old_profile <- options(rMVPA.profile_searchlight = TRUE)
  on.exit(options(old_profile), add = TRUE)

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  plan_started <- proc.time()[["elapsed"]]
  configured_backend <- configure_future_plan(config$future_backend, config$workers)
  plan_setup_seconds <- proc.time()[["elapsed"]] - plan_started

  fixture_started <- proc.time()[["elapsed"]]
  fixture <- build_mvpa_spec(config)
  fixture_setup_seconds <- proc.time()[["elapsed"]] - fixture_started
  batch_size <- if (identical(config$batch_size, "auto")) NULL else as.integer(config$batch_size)

  analysis_started <- proc.time()[["elapsed"]]
  if (identical(config$analysis, "regional")) {
    region_mask <- build_region_mask(fixture$dataset, n_regions = config$n_regions)
    args <- list(
      model_spec = fixture$model_spec,
      region_mask = region_mask,
      verbose = FALSE,
      backend = config$data_backend,
      preflight = "off"
    )
    if (!is.null(batch_size)) {
      args$batch_size <- batch_size
    }
    res <- do.call(run_regional, args)
    perf_rows <- nrow(res$performance_table %||% tibble::tibble())
    prediction_rows <- nrow(res$prediction_table %||% tibble::tibble())
    fit_count <- length(res$fits %||% list())
  } else if (identical(config$analysis, "searchlight")) {
    args <- list(
      model_spec = fixture$model_spec,
      radius = config$radius,
      method = "standard",
      verbose = FALSE,
      backend = config$data_backend,
      preflight = "off"
    )
    if (!is.null(batch_size)) {
      args$batch_size <- batch_size
    }
    res <- do.call(run_searchlight, args)
    perf_rows <- length(res$results %||% list())
    prediction_rows <- NA_integer_
    fit_count <- NA_integer_
  } else {
    stop(sprintf("Unsupported analysis '%s'", config$analysis), call. = FALSE)
  }

  analysis_seconds <- proc.time()[["elapsed"]] - analysis_started
  canonical <- canonical_result(res, config$analysis)
  if (length(canonical$values) < 1L || any(!is.finite(canonical$values))) {
    stop("Probe produced no finite canonical result values.", call. = FALSE)
  }
  receipt <- profile_receipt(res)
  signature_payload <- list(keys = canonical$keys, values = canonical$values)

  list(
    result_class = paste(class(res), collapse = ";"),
    perf_rows = perf_rows,
    prediction_rows = prediction_rows,
    fit_count = fit_count,
    configured_backend = configured_backend,
    configured_workers = as.integer(future::nbrOfWorkers()),
    analysis_seconds = as.numeric(analysis_seconds),
    plan_setup_seconds = as.numeric(plan_setup_seconds),
    fixture_setup_seconds = as.numeric(fixture_setup_seconds),
    result_n = length(canonical$values),
    result_keys = paste(canonical$keys, collapse = ";"),
    result_values = paste(
      formatC(canonical$values, digits = 17L, format = "g"),
      collapse = ";"
    ),
    result_signature = md5_object(signature_payload),
    data_bytes = as.numeric(utils::object.size(fixture$model_spec$dataset)),
    frame_bytes_sum = receipt$frame_bytes_sum,
    frame_bytes_max = receipt$frame_bytes_max,
    roi_extract_seconds = receipt$roi_extract_seconds,
    run_future_seconds = receipt$run_future_seconds,
    n_batches = receipt$n_batches,
    package_version = as.character(utils::packageVersion("rMVPA")),
    package_path = normalizePath(system.file(package = "rMVPA"), winslash = "/", mustWork = FALSE),
    future_version = as.character(utils::packageVersion("future")),
    furrr_version = as.character(utils::packageVersion("furrr")),
    shard_version = if (requireNamespace("shard", quietly = TRUE)) {
      as.character(utils::packageVersion("shard"))
    } else {
      NA_character_
    }
  )
}

config_from_env <- function() {
  list(
    analysis = env_choice("RMVPA_HPC_SWEEP_ANALYSIS", c("regional", "searchlight"), "regional"),
    future_backend = env_choice(
      "RMVPA_HPC_SWEEP_FUTURE_BACKEND",
      c("sequential", "multisession", "multicore", "callr", "mirai_multisession"),
      "sequential"
    ),
    data_backend = env_choice(
      "RMVPA_HPC_SWEEP_DATA_BACKEND",
      c("default", "shard"),
      "default"
    ),
    load_mode = env_choice(
      "RMVPA_HPC_SWEEP_LOAD_MODE",
      c("source", "installed"),
      "source"
    ),
    workers = env_int("RMVPA_HPC_SWEEP_WORKERS", 1L),
    omp_threads = env_int("RMVPA_HPC_SWEEP_OMP_THREADS", 1L),
    blas_threads = env_int("RMVPA_HPC_SWEEP_BLAS_THREADS", 1L),
    batch_size = Sys.getenv("RMVPA_HPC_SWEEP_BATCH_SIZE", "auto"),
    repetition = env_int("RMVPA_HPC_SWEEP_REPETITION", 1L),
    timeout_seconds = env_int("RMVPA_HPC_SWEEP_TIMEOUT_SECONDS", 120L),
    radius = env_int("RMVPA_HPC_SWEEP_RADIUS", 3L),
    dims = parse_dims(Sys.getenv("RMVPA_HPC_SWEEP_DIMS", "7,7,7")),
    nobs = env_int("RMVPA_HPC_SWEEP_NOBS", 48L),
    nlevels = env_int("RMVPA_HPC_SWEEP_NLEVELS", 3L),
    blocks = env_int("RMVPA_HPC_SWEEP_BLOCKS", 4L),
    n_regions = env_int("RMVPA_HPC_SWEEP_N_REGIONS", 18L),
    model = Sys.getenv("RMVPA_HPC_SWEEP_MODEL", "sda_notune"),
    log_file = Sys.getenv("RMVPA_HPC_SWEEP_LOG_FILE", ""),
    result_file = Sys.getenv("RMVPA_HPC_SWEEP_RESULT_FILE", "")
  )
}

config_id <- function(config) {
  safe_slug(sprintf(
    "%s-%s-%s-w%d-omp%d-blas%d-batch%s",
    config$analysis,
    config$future_backend,
    config$data_backend,
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
  result <- tryCatch({
    set.seed(9000L + config$repetition)
    run_probe_case(config)
  }, error = function(e) {
    status <<- "error"
    message_text <<- conditionMessage(e)
    NULL
  })

  elapsed <- as.numeric(difftime(Sys.time(), started_at, units = "secs"))
  result_value <- function(name, default) {
    if (is.null(result)) default else result[[name]] %||% default
  }

  row <- data.frame(
    config_id = config_id(config),
    run_id = run_id(config),
    analysis = config$analysis,
    future_backend = config$future_backend,
    configured_backend = as.character(result_value("configured_backend", NA_character_)),
    data_backend = config$data_backend,
    load_mode = config$load_mode,
    workers = as.integer(config$workers),
    configured_workers = as.integer(result_value("configured_workers", NA_integer_)),
    omp_threads = as.integer(config$omp_threads),
    blas_threads = as.integer(config$blas_threads),
    batch_size = as.character(config$batch_size),
    repetition = as.integer(config$repetition),
    timeout_seconds = as.integer(config$timeout_seconds),
    dims = paste(config$dims, collapse = "x"),
    nobs = as.integer(config$nobs),
    nlevels = as.integer(config$nlevels),
    blocks = as.integer(config$blocks),
    n_regions = as.integer(config$n_regions),
    model = config$model,
    status = status,
    elapsed_seconds = elapsed,
    analysis_seconds = as.numeric(result_value("analysis_seconds", NA_real_)),
    plan_setup_seconds = as.numeric(result_value("plan_setup_seconds", NA_real_)),
    fixture_setup_seconds = as.numeric(result_value("fixture_setup_seconds", NA_real_)),
    perf_rows = as.integer(result_value("perf_rows", NA_integer_)),
    prediction_rows = as.integer(result_value("prediction_rows", NA_integer_)),
    fit_count = as.integer(result_value("fit_count", NA_integer_)),
    result_class = as.character(result_value("result_class", NA_character_)),
    result_n = as.integer(result_value("result_n", NA_integer_)),
    result_keys = as.character(result_value("result_keys", NA_character_)),
    result_values = as.character(result_value("result_values", NA_character_)),
    result_signature = as.character(result_value("result_signature", NA_character_)),
    data_bytes = as.numeric(result_value("data_bytes", NA_real_)),
    frame_bytes_sum = as.numeric(result_value("frame_bytes_sum", NA_real_)),
    frame_bytes_max = as.numeric(result_value("frame_bytes_max", NA_real_)),
    roi_extract_seconds = as.numeric(result_value("roi_extract_seconds", NA_real_)),
    run_future_seconds = as.numeric(result_value("run_future_seconds", NA_real_)),
    n_batches = as.integer(result_value("n_batches", NA_integer_)),
    peak_tree_rss_bytes = NA_real_,
    peak_process_count = NA_integer_,
    memory_samples = NA_integer_,
    package_version = as.character(result_value("package_version", NA_character_)),
    package_path = as.character(result_value("package_path", NA_character_)),
    future_version = as.character(result_value("future_version", NA_character_)),
    furrr_version = as.character(result_value("furrr_version", NA_character_)),
    shard_version = as.character(result_value("shard_version", NA_character_)),
    r_version = R.version.string,
    platform = R.version$platform,
    available_cores = as.integer(tryCatch(
      parallelly::availableCores(),
      error = function(e) NA_integer_
    )),
    message = message_text,
    pid = as.integer(pid),
    host = host,
    started_at = timestamp_iso(started_at),
    finished_at = timestamp_iso(Sys.time()),
    log_file = config$log_file,
    stringsAsFactors = FALSE
  )

  if (nzchar(config$result_file)) {
    atomic_write_csv(row, config$result_file)
  } else {
    print(row)
  }
  invisible(row)
}

expand_configs <- function() {
  analyses <- env_csv("RMVPA_HPC_SWEEP_ANALYSES", default = "regional")
  backends <- env_csv("RMVPA_HPC_SWEEP_BACKENDS", default = c("sequential", "multisession", "multicore"))
  data_backends <- env_csv("RMVPA_HPC_SWEEP_DATA_BACKENDS", default = "default")
  workers <- parse_positive_ints(env_csv("RMVPA_HPC_SWEEP_WORKER_COUNTS", default = c("1", "2")), "RMVPA_HPC_SWEEP_WORKER_COUNTS")
  omp_threads <- parse_positive_ints(env_csv("RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS", default = c("1", "2")), "RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS")
  blas_threads <- parse_positive_ints(env_csv("RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS", default = c("1")), "RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS")
  batch_sizes <- parse_batch_values(env_csv("RMVPA_HPC_SWEEP_BATCH_SIZES", default = "auto"))
  repetitions <- seq_len(env_int("RMVPA_HPC_SWEEP_REP", 1L))
  timeout_seconds <- env_int("RMVPA_HPC_SWEEP_TIMEOUT_SECONDS", 120L)
  radius <- env_int("RMVPA_HPC_SWEEP_RADIUS", 3L)
  load_mode <- env_choice(
    "RMVPA_HPC_SWEEP_LOAD_MODE",
    c("source", "installed"),
    "source"
  )
  dims <- parse_dims(Sys.getenv("RMVPA_HPC_SWEEP_DIMS", "7,7,7"))
  nobs <- env_int("RMVPA_HPC_SWEEP_NOBS", 48L)
  nlevels <- env_int("RMVPA_HPC_SWEEP_NLEVELS", 3L)
  blocks <- env_int("RMVPA_HPC_SWEEP_BLOCKS", 4L)
  n_regions <- env_int("RMVPA_HPC_SWEEP_N_REGIONS", 18L)
  model <- Sys.getenv("RMVPA_HPC_SWEEP_MODEL", "sda_notune")

  grid <- expand.grid(
    analysis = analyses,
    future_backend = backends,
    data_backend = data_backends,
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
  grid$load_mode <- load_mode
  grid$dims <- paste(dims, collapse = ",")
  grid$nobs <- nobs
  grid$nlevels <- nlevels
  grid$blocks <- blocks
  grid$n_regions <- n_regions
  grid$model <- model
  grid
}

launch_child <- function(script_path, config, log_file, result_file) {
  env <- c(
    "RMVPA_HPC_SWEEP_MODE=child",
    sprintf("RMVPA_HPC_SWEEP_ANALYSIS=%s", config$analysis),
    sprintf("RMVPA_HPC_SWEEP_FUTURE_BACKEND=%s", config$future_backend),
    sprintf("RMVPA_HPC_SWEEP_DATA_BACKEND=%s", config$data_backend),
    sprintf("RMVPA_HPC_SWEEP_LOAD_MODE=%s", config$load_mode),
    sprintf("RMVPA_HPC_SWEEP_WORKERS=%d", as.integer(config$workers)),
    sprintf("RMVPA_HPC_SWEEP_OMP_THREADS=%d", as.integer(config$omp_threads)),
    sprintf("RMVPA_HPC_SWEEP_BLAS_THREADS=%d", as.integer(config$blas_threads)),
    sprintf("RMVPA_HPC_SWEEP_BATCH_SIZE=%s", config$batch_size),
    sprintf("RMVPA_HPC_SWEEP_REPETITION=%d", as.integer(config$repetition)),
    sprintf("RMVPA_HPC_SWEEP_TIMEOUT_SECONDS=%d", as.integer(config$timeout_seconds)),
    sprintf("RMVPA_HPC_SWEEP_RADIUS=%d", as.integer(config$radius)),
    sprintf("RMVPA_HPC_SWEEP_DIMS=%s", paste(config$dims, collapse = ",")),
    sprintf("RMVPA_HPC_SWEEP_NOBS=%d", as.integer(config$nobs)),
    sprintf("RMVPA_HPC_SWEEP_NLEVELS=%d", as.integer(config$nlevels)),
    sprintf("RMVPA_HPC_SWEEP_BLOCKS=%d", as.integer(config$blocks)),
    sprintf("RMVPA_HPC_SWEEP_N_REGIONS=%d", as.integer(config$n_regions)),
    sprintf("RMVPA_HPC_SWEEP_MODEL=%s", config$model),
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

process_tree_sample <- function(pid) {
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(list(rss_bytes = NA_real_, process_count = NA_integer_))
  }
  root <- tryCatch(ps::ps_handle(as.integer(pid)), error = function(e) NULL)
  if (is.null(root)) {
    list(rss_bytes = NA_real_, process_count = NA_integer_)
  } else {
    children <- tryCatch(
      ps::ps_children(root, recursive = TRUE),
      error = function(e) list()
    )
    handles <- c(list(root), children)
    rss <- vapply(handles, function(handle) {
      tryCatch(
        as.numeric(ps::ps_memory_info(handle)[["rss"]]),
        error = function(e) NA_real_
      )
    }, numeric(1))
    keep <- is.finite(rss)
    list(
      rss_bytes = if (any(keep)) sum(rss[keep]) else NA_real_,
      process_count = if (any(keep)) as.integer(sum(keep)) else NA_integer_
    )
  }
}

wait_for_result <- function(pid, result_file, timeout_seconds) {
  started <- Sys.time()
  peak_tree_rss_bytes <- NA_real_
  peak_process_count <- NA_integer_
  memory_samples <- 0L
  repeat {
    sample <- process_tree_sample(pid)
    if (is.finite(sample$rss_bytes)) {
      peak_tree_rss_bytes <- if (is.finite(peak_tree_rss_bytes)) {
        max(peak_tree_rss_bytes, sample$rss_bytes)
      } else {
        sample$rss_bytes
      }
      peak_process_count <- if (is.na(peak_process_count)) {
        sample$process_count
      } else {
        max(peak_process_count, sample$process_count, na.rm = TRUE)
      }
      memory_samples <- memory_samples + 1L
    }
    if (result_file_ready(result_file)) {
      state <- "finished"
      break
    }
    if (!is_process_running(pid)) {
      state <- "exited"
      break
    }
    if (as.numeric(difftime(Sys.time(), started, units = "secs")) > timeout_seconds) {
      state <- "timeout"
      break
    }
    Sys.sleep(0.25)
  }
  list(
    state = state,
    peak_tree_rss_bytes = peak_tree_rss_bytes,
    peak_process_count = peak_process_count,
    memory_samples = memory_samples
  )
}

row_from_config <- function(config, status, log_file, message = NA_character_) {
  data.frame(
    config_id = config_id(config),
    run_id = run_id(config),
    analysis = config$analysis,
    future_backend = config$future_backend,
    configured_backend = NA_character_,
    data_backend = config$data_backend,
    load_mode = config$load_mode,
    workers = as.integer(config$workers),
    configured_workers = NA_integer_,
    omp_threads = as.integer(config$omp_threads),
    blas_threads = as.integer(config$blas_threads),
    batch_size = as.character(config$batch_size),
    repetition = as.integer(config$repetition),
    timeout_seconds = as.integer(config$timeout_seconds),
    dims = gsub(",", "x", paste(config$dims, collapse = ","), fixed = TRUE),
    nobs = as.integer(config$nobs),
    nlevels = as.integer(config$nlevels),
    blocks = as.integer(config$blocks),
    n_regions = as.integer(config$n_regions),
    model = config$model,
    status = status,
    elapsed_seconds = NA_real_,
    analysis_seconds = NA_real_,
    plan_setup_seconds = NA_real_,
    fixture_setup_seconds = NA_real_,
    perf_rows = NA_integer_,
    prediction_rows = NA_integer_,
    fit_count = NA_integer_,
    result_class = NA_character_,
    result_n = NA_integer_,
    result_keys = NA_character_,
    result_values = NA_character_,
    result_signature = NA_character_,
    data_bytes = NA_real_,
    frame_bytes_sum = NA_real_,
    frame_bytes_max = NA_real_,
    roi_extract_seconds = NA_real_,
    run_future_seconds = NA_real_,
    n_batches = NA_integer_,
    peak_tree_rss_bytes = NA_real_,
    peak_process_count = NA_integer_,
    memory_samples = NA_integer_,
    package_version = NA_character_,
    package_path = NA_character_,
    future_version = NA_character_,
    furrr_version = NA_character_,
    shard_version = NA_character_,
    r_version = R.version.string,
    platform = R.version$platform,
    available_cores = as.integer(tryCatch(
      parallelly::availableCores(),
      error = function(e) NA_integer_
    )),
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
    median_ok <- function(x) {
      keep <- ok & is.finite(x)
      if (any(keep)) stats::median(x[keep]) else NA_real_
    }
    problem_idx <- which(!is.na(df$message) & nzchar(df$message))
    first_problem <- if (length(problem_idx) > 0L) df$message[[problem_idx[[1L]]]] else NA_character_
    data.frame(
      config_id = df$config_id[[1L]],
      analysis = df$analysis[[1L]],
      future_backend = df$future_backend[[1L]],
      data_backend = df$data_backend[[1L]],
      load_mode = df$load_mode[[1L]],
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
      median_elapsed_ok = median_ok(df$elapsed_seconds),
      median_analysis_seconds_ok = median_ok(df$analysis_seconds),
      median_plan_setup_seconds_ok = median_ok(df$plan_setup_seconds),
      median_frame_bytes_max_ok = median_ok(df$frame_bytes_max),
      median_roi_extract_seconds_ok = median_ok(df$roi_extract_seconds),
      median_run_future_seconds_ok = median_ok(df$run_future_seconds),
      median_peak_tree_rss_bytes_ok = median_ok(df$peak_tree_rss_bytes),
      example_message = first_problem,
      example_log_file = df$log_file[[1L]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

config_availability <- function(config) {
  future_status <- backend_availability(config$future_backend)
  if (!isTRUE(future_status$ok)) return(future_status)
  data_backend_availability(config$data_backend)
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
      avail <- config_availability(cfg)
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
    avail <- config_availability(cfg)

    if (!isTRUE(avail$ok)) {
      rows[[i]] <- row_from_config(cfg, status = "skip", log_file = log_file, message = avail$reason)
      next
    }

    unlink(result_file, force = TRUE)
    message(sprintf(
      "[%d/%d] %s future=%s data=%s workers=%d omp=%d blas=%d batch=%s",
      i, nrow(configs), cfg$analysis, cfg$future_backend, cfg$data_backend, cfg$workers,
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

    wait <- wait_for_result(pid, result_file, timeout_seconds = cfg$timeout_seconds)

    if (identical(wait$state, "timeout")) {
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
      row$peak_tree_rss_bytes <- wait$peak_tree_rss_bytes
      row$peak_process_count <- wait$peak_process_count
      row$memory_samples <- wait$memory_samples
      rows[[i]] <- row
      next
    }

    if (file.exists(result_file)) {
      row <- read_result_row(result_file)
      if (inherits(row, "error")) {
        msg <- sprintf(
          "Child published an invalid result row (%s). Log tail:\n%s",
          conditionMessage(row),
          tail_log(log_file, n = 20L)
        )
        row <- row_from_config(
          cfg, status = "crash", log_file = log_file, message = msg
        )
        row$pid <- as.integer(pid)
        row$started_at <- timestamp_iso(started)
        row$finished_at <- timestamp_iso(Sys.time())
      } else {
        row$log_file <- log_file
      }
      row$peak_tree_rss_bytes <- wait$peak_tree_rss_bytes
      row$peak_process_count <- wait$peak_process_count
      row$memory_samples <- wait$memory_samples
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
    row$peak_tree_rss_bytes <- wait$peak_tree_rss_bytes
    row$peak_process_count <- wait$peak_process_count
    row$memory_samples <- wait$memory_samples
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

if (sys.nframe() == 0L) {
  main()
}
