#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
})

capability_patterns <- list(
  core = paste(
    c(
      "test_validate_analysis",
      "test_output_schema",
      "test_save_results",
      "test_fit_roi",
      "test_plugin_helpers",
      "test_plugin_extension_api",
      "test_searchlight_engine_audit",
      "test_dataset",
      "test_average_labels",
      "test-crossv_functions",
      "test_balance_partition",
      "test_distfun",
      "test_table_to_rdm"
    ),
    collapse = "|"
  ),
  searchlight = paste(
    c(
      "test_mvpa_searchlight",
      "test_searchlight_",
      "test_clustered_searchlight",
      "test_swift_searchlight",
      "test_dual_lda_fast_searchlight",
      "test_filter_roi_fast_backend"
    ),
    collapse = "|"
  ),
  regional = paste(
    c(
      "test_mvpa_regional",
      "test_rsa_regional",
      "test_vector_rsa_regional",
      "test_global_analysis",
      "test_region_importance",
      "test_banded_ridge_da_model"
    ),
    collapse = "|"
  ),
  shard = paste(
    c(
      "test_feature_rsa_shard",
      "test_shard_backend"
    ),
    collapse = "|"
  ),
  remap = paste(
    c(
      "test_remap_",
      "test-repmap_model",
      "test-repmed_model",
      "test-repnet_model",
      "test_naive_xdec_model",
      "test_subspace_alignment_model"
    ),
    collapse = "|"
  )
)

args <- commandArgs(trailingOnly = TRUE)
capability <- if (length(args) == 0L) "core" else args[[1L]]

if (!capability %in% names(capability_patterns)) {
  stop(
    sprintf(
      "Unknown capability '%s'. Expected one of: %s",
      capability,
      paste(names(capability_patterns), collapse = ", ")
    ),
    call. = FALSE
  )
}

pattern <- capability_patterns[[capability]]
test_files <- list.files("tests/testthat", pattern = "^test.*\\.R$", full.names = FALSE)
matched <- grep(pattern, test_files, value = TRUE)

if (length(matched) == 0L) {
  stop("Capability filter matched zero test files: ", capability, call. = FALSE)
}

message(sprintf("Running capability '%s' with %d test file(s).", capability, length(matched)))
message(paste("Matched files:", paste(sort(matched), collapse = ", ")))

invisible(lapply(
  file.path("tests", "testthat", matched),
  function(path) {
    testthat::test_file(
      path,
      stop_on_failure = TRUE,
      load_package = "source"
    )
  }
))
