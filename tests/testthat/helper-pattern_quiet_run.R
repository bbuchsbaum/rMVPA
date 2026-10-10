# Evaluate an expression with messages suppressed and stdout captured (pattern
# model tests). helper-test-logging.R already pins the log threshold; this only
# hides the iteration banners that the runners print to stdout.

quiet_run <- function(expr) {
  utils::capture.output(res <- suppressMessages(expr))
  res
}
