#' List the rMVPA API Lifecycle Registry
#'
#' Returns a compact registry describing which exported entry points are
#' considered \code{"stable"}, \code{"experimental"}, or
#' \code{"developer"}.
#'
#' @return A data frame with columns \code{symbol}, \code{lifecycle}, and
#'   \code{notes}.
#' @export
rmvpa_api_lifecycle <- function() {
  data.frame(
    symbol = c(
      "mvpa_config",
      "build_analysis",
      "run_analysis",
      "save_results",
      "install_cli",
      "validate_analysis",
      "mvpa_dataset",
      "mvpa_design",
      "load_model",
      "mvpa_model",
      "run_searchlight",
      "run_regional",
      "run_global",
      "searchlight_engines",
      "explain_searchlight_engine",
      "validate_plugin_model",
      "validate_model_spec",
      "mock_roi_data",
      "mock_context",
      "cv_evaluate_roi",
      "use_shard"
    ),
    lifecycle = c(
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "stable",
      "experimental",
      "experimental",
      "developer",
      "developer",
      "developer",
      "developer",
      "developer",
      "experimental"
    ),
    notes = c(
      "Public workflow configuration constructor.",
      "Public analysis builder for CLI and scripted use.",
      "Public workflow runner with optional automatic preflight checks.",
      "Stable result writer and reproducibility bundle entry point.",
      "Stable installer for packaged CLI wrappers.",
      "Stable static methodology validator.",
      "Stable dataset constructor.",
      "Stable design constructor.",
      "Stable model-registry lookup.",
      "Stable core model-spec constructor.",
      "Stable searchlight runner.",
      "Stable regional runner.",
      "Stable global runner.",
      "Fast-path engine registry summary.",
      "Eligibility and selection report for searchlight engines.",
      "Plugin contract checker for one ROI.",
      "Developer lint for plugin model specs.",
      "Developer helper for plugin/unit tests.",
      "Developer helper for plugin/unit tests.",
      "Developer wrapper around ROI CV helpers.",
      "Experimental shared-memory backend selector."
    ),
    stringsAsFactors = FALSE
  )
}

#' Return the Stable rMVPA API Surface
#'
#' @return A character vector of exported symbols designated as stable.
#' @export
rmvpa_stable_api <- function() {
  reg <- rmvpa_api_lifecycle()
  reg$symbol[reg$lifecycle == "stable"]
}
