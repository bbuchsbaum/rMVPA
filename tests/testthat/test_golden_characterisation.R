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
      actual <- golden_run_scenario(nm)
      expected <- fixture$value
      if (identical(nm, "regional_sda_notune")) {
        # SDA rounds its posterior before wrap_result() normalizes each row.
        # A different summation precision can change the normalized values by
        # 2e-16, breaking exact ties and moving AUC by a discrete pair count.
        # Keep frozen probabilities and every other field at the same tolerance;
        # compare the full result to the original estimator on this platform.
        reference <- golden_external_sda()
        expect_equal(actual, reference, tolerance = golden_tolerance,
                     info = "native SDA versus original sda::sda")
        expect_equal(actual$performance$AUC, golden_pairwise_auc(actual),
                     tolerance = golden_tolerance, info = "native SDA pair-count AUC")
        expect_equal(expected$performance$AUC, golden_pairwise_auc(expected),
                     tolerance = golden_tolerance, info = "frozen SDA pair-count AUC")
        actual$performance$AUC <- NULL
        expected$performance$AUC <- NULL
      }
      expect_equal(actual, expected,
                   tolerance = golden_tolerance, info = nm)
    })
  })
}
