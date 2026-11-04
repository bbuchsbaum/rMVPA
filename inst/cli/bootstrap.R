rmvpa_cli_boot <- function(mode, package_root = NULL) {
  # mode "dsl" is the `rmvpa` command; the others are the legacy wrappers.
  call_cli <- function() {
    ns <- asNamespace("rMVPA")
    if (identical(mode, "dsl")) {
      return(get("rmvpa_cli", envir = ns, inherits = FALSE)(commandArgs(trailingOnly = TRUE)))
    }
    get("run_cli_command", envir = ns, inherits = FALSE)(mode)
    invisible(0L)
  }
  if (!is.null(package_root) &&
      file.exists(file.path(package_root, "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    suppressWarnings(
      suppressPackageStartupMessages(
        pkgload::load_all(package_root, export_all = FALSE, helpers = FALSE, quiet = TRUE)
      )
    )
    return(call_cli())
  }

  if (requireNamespace("rMVPA", quietly = TRUE)) {
    return(call_cli())
  }

  stop(
    "Could not load rMVPA. Install the package first, or install 'pkgload' to run the CLI from a source checkout.",
    call. = FALSE
  )
}
