suppressPackageStartupMessages({
  library(testthat)
})

env_truthy <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) {
    return(isTRUE(default))
  }
  tolower(raw) %in% c("1", "true", "t", "yes", "y", "on")
}

set_default_env <- function(name, value) {
  if (!nzchar(Sys.getenv(name, ""))) {
    do.call(Sys.setenv, setNames(list(value), name))
  }
}

is_process_running <- function(pid) {
  pid <- suppressWarnings(as.integer(pid))
  if (!is.finite(pid) || is.na(pid) || pid <= 0L) {
    return(FALSE)
  }
  if (.Platform$OS.type == "windows") {
    # Cross-platform stale-lock handling on Windows is harder without extra deps.
    # Keep lock conservative and rely on stale timeout there.
    return(TRUE)
  }
  status <- suppressWarnings(
    system2("kill", c("-0", as.character(pid)), stdout = FALSE, stderr = FALSE)
  )
  identical(as.integer(status), 0L)
}

acquire_perf_guardrail_lock <- function(lock_dir, stale_seconds = 6 * 3600) {
  lock_dir <- normalizePath(lock_dir, winslash = "/", mustWork = FALSE)
  lock_parent <- dirname(lock_dir)
  dir.create(lock_parent, recursive = TRUE, showWarnings = FALSE)

  try_make <- function() dir.create(lock_dir, showWarnings = FALSE)
  if (isTRUE(try_make())) {
    return(lock_dir)
  }

  owner_file <- file.path(lock_dir, "owner")
  lock_pid <- NA_integer_
  lock_started <- NA_character_
  if (file.exists(owner_file)) {
    lock_lines <- readLines(owner_file, warn = FALSE)
    if (length(lock_lines) >= 1L) {
      lock_pid <- suppressWarnings(as.integer(lock_lines[[1L]]))
    }
    if (length(lock_lines) >= 2L) {
      lock_started <- lock_lines[[2L]]
    }
  }

  lock_age <- Inf
  info <- suppressWarnings(file.info(lock_dir))
  if (!is.null(info) && nrow(info) == 1L && !is.na(info$mtime)) {
    lock_age <- as.numeric(difftime(Sys.time(), info$mtime, units = "secs"))
  }

  stale_by_pid <- !is.na(lock_pid) && !is_process_running(lock_pid)
  stale_by_age <- is.finite(lock_age) && (lock_age > stale_seconds)
  if (isTRUE(stale_by_pid) || isTRUE(stale_by_age)) {
    message(
      sprintf(
        "Removing stale perf lock at %s (pid=%s started=%s age=%.0fs).",
        lock_dir,
        ifelse(is.na(lock_pid), "unknown", as.character(lock_pid)),
        ifelse(is.na(lock_started), "unknown", lock_started),
        ifelse(is.finite(lock_age), lock_age, -1)
      )
    )
    unlink(lock_dir, recursive = TRUE, force = TRUE)
    if (isTRUE(try_make())) {
      return(lock_dir)
    }
  }

  stop(
    sprintf(
      paste0(
        "Perf guardrails lock is already held: %s. ",
        "Set RMVPA_PERF_SINGLETON=false to bypass lock if you intentionally want parallel runs."
      ),
      lock_dir
    ),
    call. = FALSE
  )
}

release_perf_guardrail_lock <- function(lock_dir) {
  unlink(lock_dir, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

set_default_env("RGL_USE_NULL", "TRUE")
set_default_env("NOT_CRAN", "true")
set_default_env("RMVPA_RUN_PERF_TESTS", "true")
set_default_env("RMVPA_PERF_REP", "1")
set_default_env("RMVPA_RUN_NIGHTLY_PERF_TESTS", "false")
set_default_env("RMVPA_NIGHTLY_PERF_REP", "2")
set_default_env("R_PROGRESSR_ENABLE", "false")
set_default_env("TESTTHAT_CPUS", "1")

main <- function() {
  default_filters <- c(
    "searchlight_framework_perf",
    "clustered_searchlight_perf",
    "dual_lda_fast_searchlight_perf",
    "swift_searchlight_perf",
    "rsa_fast_kernel_perf",
    "naive_xdec_fast_kernel_perf"
  )
  perf_filter <- Sys.getenv("RMVPA_PERF_FILTER", "")
  if (!nzchar(perf_filter)) {
    perf_filter <- paste(default_filters, collapse = "|")
  }

  if (env_truthy("RMVPA_PERF_SINGLETON", default = TRUE)) {
    lock_dir <- Sys.getenv(
      "RMVPA_PERF_LOCK_DIR",
      file.path(getwd(), ".tmp", "perf_guardrails.lock")
    )
    stale_seconds <- suppressWarnings(as.numeric(Sys.getenv("RMVPA_PERF_LOCK_STALE_SECONDS", "21600")))
    if (!is.finite(stale_seconds) || is.na(stale_seconds) || stale_seconds <= 0) {
      stale_seconds <- 21600
    }
    lock_dir <- acquire_perf_guardrail_lock(lock_dir, stale_seconds = stale_seconds)
    writeLines(
      c(
        as.character(Sys.getpid()),
        format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      ),
      con = file.path(lock_dir, "owner")
    )
    on.exit(release_perf_guardrail_lock(lock_dir), add = TRUE)
  }

  message(
    sprintf(
      paste0(
        "Running perf guardrails with filter=%s ",
        "RGL_USE_NULL=%s NOT_CRAN=%s RMVPA_RUN_PERF_TESTS=%s RMVPA_PERF_REP=%s ",
        "RMVPA_RUN_NIGHTLY_PERF_TESTS=%s RMVPA_NIGHTLY_PERF_REP=%s TESTTHAT_CPUS=%s"
      ),
      perf_filter,
      Sys.getenv("RGL_USE_NULL", ""),
      Sys.getenv("NOT_CRAN", ""),
      Sys.getenv("RMVPA_RUN_PERF_TESTS", ""),
      Sys.getenv("RMVPA_PERF_REP", ""),
      Sys.getenv("RMVPA_RUN_NIGHTLY_PERF_TESTS", ""),
      Sys.getenv("RMVPA_NIGHTLY_PERF_REP", ""),
      Sys.getenv("TESTTHAT_CPUS", "")
    )
  )

  testthat::test_local(
    ".",
    filter = perf_filter,
    stop_on_failure = TRUE,
    load_package = "source"
  )
}

main()
