# Helper file to control logging during tests
# This file is automatically loaded before tests run

# Suppress INFO and WARN level logging during tests
# Only show ERROR and FATAL messages
if (requireNamespace("futile.logger", quietly = TRUE)) {
  # Set threshold for root logger and rMVPA namespace
  futile.logger::flog.threshold(futile.logger::ERROR, name = "ROOT")
  futile.logger::flog.threshold(futile.logger::ERROR, name = "rMVPA")
  
  # Also set for the default logger
  futile.logger::flog.threshold(futile.logger::ERROR)
}

# Optionally allow override via environment variable for debugging
log_level_env <- Sys.getenv("RMVPA_TEST_LOG_LEVEL", "")
if (log_level_env != "") {
  log_level <- switch(log_level_env,
    "DEBUG" = futile.logger::DEBUG,
    "INFO" = futile.logger::INFO,
    "WARN" = futile.logger::WARN,
    "ERROR" = futile.logger::ERROR,
    "FATAL" = futile.logger::FATAL,
    futile.logger::ERROR  # default
  )
  futile.logger::flog.threshold(log_level, name = "ROOT")
  futile.logger::flog.threshold(log_level, name = "rMVPA")
  futile.logger::flog.threshold(log_level)
}

# Suppress package version warnings during library() calls
# These are harmless but create visual noise
options(warn = -1)  # Temporarily suppress warnings
options(progress_enabled = FALSE)
options(cli.progress_show_after = Inf)
Sys.setenv(R_PROGRESSR_ENABLE = "false")

# Load commonly used packages with suppressed warnings
suppressWarnings({
  library(neuroim2)
  library(neurosurf)
  library(testthat)
  library(assertthat)
})

# Restore normal warning behavior
options(warn = 0)
