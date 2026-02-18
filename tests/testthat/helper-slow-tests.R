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
