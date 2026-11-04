rmvpa_cli_boot <- function(mode, package_root = NULL) {
  if (!is.null(package_root) &&
      file.exists(file.path(package_root, "DESCRIPTION")) &&
      requireNamespace("pkgload", quietly = TRUE)) {
    suppressWarnings(
      suppressPackageStartupMessages(
        pkgload::load_all(package_root, export_all = FALSE, helpers = FALSE, quiet = TRUE)
      )
    )
    get("run_cli_command", envir = asNamespace("rMVPA"), inherits = FALSE)(mode)
    return(invisible(NULL))
  }

  if (requireNamespace("rMVPA", quietly = TRUE)) {
    get("run_cli_command", envir = asNamespace("rMVPA"), inherits = FALSE)(mode)
    return(invisible(NULL))
  }

  stop(
    "Could not load rMVPA. Install the package first, or install 'pkgload' to run the CLI from a source checkout.",
    call. = FALSE
  )
}
