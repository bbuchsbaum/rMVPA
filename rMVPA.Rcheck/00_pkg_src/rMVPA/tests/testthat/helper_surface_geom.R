rmvpa_test_surface_geom_file <- function() {
  # Prefer neurosurf's built-in small geometry used elsewhere in the package.
  if (requireNamespace("neurosurf", quietly = TRUE)) {
    fname <- system.file("extdata/std.8_lh.inflated.asc", package = "neurosurf")
    if (!identical(fname, "") && file.exists(fname)) {
      return(fname)
    }
  }

  # Back-compat fallback for older test environments.
  if (requireNamespace("neuroim2", quietly = TRUE)) {
    fname <- system.file("extdata/std.lh.smoothwm.asc", package = "neuroim2")
    if (!identical(fname, "") && file.exists(fname)) {
      return(fname)
    }
  }

  ""
}

