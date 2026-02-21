rmvpa_extended_tests_enabled <- function() {
  flag <- tolower(Sys.getenv("RMVPA_RUN_EXTENDED_TESTS", "false"))
  flag %in% c("1", "true", "yes", "on")
}

skip_if_not_extended_tests <- function() {
  testthat::skip_if_not(
    rmvpa_extended_tests_enabled(),
    message = "Set RMVPA_RUN_EXTENDED_TESTS=true to run long integration tests."
  )
}

rmvpa_perf_tests_enabled <- function() {
  flag <- tolower(Sys.getenv("RMVPA_RUN_PERF_TESTS", "false"))
  flag %in% c("1", "true", "yes", "on")
}

skip_if_not_perf_tests <- function() {
  testthat::skip_if_not(
    rmvpa_perf_tests_enabled(),
    message = "Set RMVPA_RUN_PERF_TESTS=true to run performance guardrail tests."
  )
}

rmvpa_nightly_perf_tests_enabled <- function() {
  flag <- tolower(Sys.getenv("RMVPA_RUN_NIGHTLY_PERF_TESTS", "false"))
  flag %in% c("1", "true", "yes", "on")
}

skip_if_not_nightly_perf_tests <- function() {
  testthat::skip_if_not(
    rmvpa_nightly_perf_tests_enabled(),
    message = "Set RMVPA_RUN_NIGHTLY_PERF_TESTS=true to run nightly performance tests."
  )
}
