#!/usr/bin/env Rscript

argv <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", argv, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
} else {
  normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = FALSE)
}

package_root <- NULL
bootstrap <- NULL

if (nzchar(script_path)) {
  root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
  candidates <- c(
    file.path(root, "inst", "cli", "bootstrap.R"),
    file.path(root, "cli", "bootstrap.R")
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) > 0L &&
      file.exists(file.path(root, "DESCRIPTION")) &&
      dir.exists(file.path(root, "scripts"))) {
    package_root <- root
    bootstrap <- hits[[1]]
  }
}

if (is.null(bootstrap)) {
  bootstrap <- system.file("cli", "bootstrap.R", package = "rMVPA")
}

if (!nzchar(bootstrap) || !file.exists(bootstrap)) {
  stop("Could not locate the rMVPA CLI bootstrap script.", call. = FALSE)
}

source(bootstrap, local = TRUE)
rmvpa_cli_boot("searchlight", package_root = package_root)
