#' rmvpadsl: a declarative analysis language for rMVPA
#'
#' Start with [rmvpa_new()] or `rmvpa new <bids_root>`, then [rmvpa_check()],
#' [rmvpa_run()] and [rmvpa_status()]. Each analysis compiles to an ordinary
#' rMVPA `model_spec`; [rmvpa_compile()] returns it for one subject.
#'
#' This package depends only on rMVPA's exported API. Facts about rMVPA's
#' behaviour that the language relies on (which scopes each model supports,
#' which permutation nulls exist) live in `.dsl_scope_support()`,
#' `.dsl_inference_support()` and `.dsl_required_packages()`, and are pinned by
#' the end-to-end tests, which run every verb against the real constructors.
#'
#' @import rMVPA
#' @keywords internal
"_PACKAGE"

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Install the `rmvpa` command
#'
#' Copies the `rmvpa` launcher into a directory on your `PATH`.
#'
#' @param dest_dir Destination directory, for example `"~/.local/bin"`.
#'   Required: nothing is written unless you choose it.
#' @param overwrite Replace an existing file.
#' @return The installed path, invisibly.
#' @export
install_rmvpa_command <- function(dest_dir, overwrite = FALSE) {
  if (missing(dest_dir)) {
    stop("install_rmvpa_command() needs `dest_dir`, a directory on your PATH (for example \"~/.local/bin\").",
         call. = FALSE)
  }
  src <- system.file("exec", "rmvpa", package = "rmvpadsl")
  if (!nzchar(src)) stop("The rmvpa launcher is missing; reinstall rmvpadsl.", call. = FALSE)
  dest_dir <- path.expand(dest_dir)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  dest <- file.path(dest_dir, "rmvpa")
  if (file.exists(dest) && !overwrite) {
    stop("Refusing to overwrite ", dest, ". Re-run with overwrite = TRUE.", call. = FALSE)
  }
  if (!file.copy(src, dest, overwrite = overwrite)) stop("Could not copy the launcher to ", dest, call. = FALSE)
  Sys.chmod(dest, mode = "755")
  if (!dest_dir %in% strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]) {
    message("Installed ", dest, ". Add ", dest_dir, " to PATH to call `rmvpa` directly.")
  }
  invisible(dest)
}
