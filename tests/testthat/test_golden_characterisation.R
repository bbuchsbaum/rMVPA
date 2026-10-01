# Golden characterisation tests: the core methods must reproduce the frozen
# outputs in fixtures/golden/. These guard refactors and performance work.
# They do not establish correctness. For a deliberate behaviour change,
# regenerate with data-raw/golden/make_golden.R and explain the change in
# NEWS.md (see helper-golden.R).

testthat::skip_if_not_installed("neuroim2")
testthat::skip_if_not_installed("sda")

golden_tolerance <- 1e-8

for (scenario in names(golden_scenarios())) {
  local({
    nm <- scenario
    test_that(sprintf("golden scenario '%s' reproduces its fixture", nm), {
      if (startsWith(nm, "sl_")) {
        # Searchlight scenarios run the general-purpose iterator over 64
        # centres; keep them out of the CRAN time budget.
        skip_on_cran()
      }
      path <- file.path(golden_fixture_dir(), paste0(nm, ".rds"))
      expect_true(file.exists(path), info = paste("missing fixture", path))
      skip_if_not(file.exists(path))
      fixture <- readRDS(path)
      expect_equal(golden_run_scenario(nm), fixture$value,
                   tolerance = golden_tolerance, info = nm)
    })
  })
}
