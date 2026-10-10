# Muffle worker-version warnings emitted by shard/feature RSA worker tests.

muffle_worker_version_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("was built under R version", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
