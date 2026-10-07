# Regenerate the golden characterisation fixtures used by
# tests/testthat/test_golden_characterisation.R.
#
# Run from the package root:
#   Rscript data-raw/golden/make_golden.R              # every scenario
#   Rscript data-raw/golden/make_golden.R name1 name2  # selected scenarios
#
# A fixture records current behaviour, not correctness. Regenerate one only
# for a deliberate behaviour change, in the same commit as that change, with a
# NEWS.md entry explaining why its values moved. Scenarios are defined in
# tests/testthat/helper-golden.R.

pkgload::load_all(".", quiet = TRUE)
library(testthat)
source(file.path("tests", "testthat", "helper-golden.R"))

out_dir <- file.path("tests", "testthat", "fixtures", "golden")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

requested <- commandArgs(trailingOnly = TRUE)
all_names <- names(golden_scenarios())
targets <- if (length(requested)) requested else all_names
unknown <- setdiff(targets, all_names)
if (length(unknown)) stop("Unknown scenario(s): ", paste(unknown, collapse = ", "))

sha <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE), error = function(e) NA_character_)
meta_common <- list(
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  git_sha = sha,
  r_version = R.version.string,
  blas = extSoftVersion()[["BLAS"]],
  neuroim2 = as.character(utils::packageVersion("neuroim2")),
  sda = tryCatch(as.character(utils::packageVersion("sda")), error = function(e) NA_character_)
)

for (nm in targets) {
  t0 <- Sys.time()
  value <- golden_run_scenario(nm)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fixture <- list(value = value, meta = c(meta_common, scenario = nm, seconds = secs))
  if (identical(nm, "rsa_spearman_pearson")) {
    fixture$distance_reference <- list(
      value = golden_spearman_rsa_components()$distances,
      meta = meta_common
    )
  }
  saveRDS(fixture,
          file.path(out_dir, paste0(nm, ".rds")), version = 2)
  cat(sprintf("%-24s %6.1fs\n", nm, secs))
}
