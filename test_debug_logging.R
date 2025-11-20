#!/usr/bin/env Rscript
# Test script to demonstrate debug logging for combine_randomized

library(rMVPA)

# Enable DEBUG logging
cat("Setting log level to DEBUG...\n")
rMVPA::set_log_level("DEBUG")

# Load the package
devtools::load_all()

cat("\n=== Running a simple randomized searchlight test ===\n\n")

# This will trigger the debug logging we just added
testthat::test_file('tests/testthat/test_mvpa_searchlight.R',
                   filter = "randomized")

cat("\n=== Test complete ===\n")
cat("\nTo use in your HPC job, add this at the start of your script:\n")
cat("  rMVPA::set_log_level('DEBUG')\n")
cat("\nThis will show detailed progress through:\n")
cat("  - Each iteration of randomized searchlight\n")
cat("  - ROI processing counts\n")
cat("  - Triplet accumulation and sparse matrix construction\n")
cat("  - Memory usage estimates\n")
