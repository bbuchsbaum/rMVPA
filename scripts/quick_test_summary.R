#!/usr/bin/env Rscript
###############################################################################
# Quick Test Summary Script
#
# Provides a quick overview of test coverage without full covr analysis
###############################################################################

cat("\n=== rMVPA Test Summary ===\n\n")

# Count R files
r_files <- list.files("R", pattern = "\\.R$", full.names = FALSE)
cat(sprintf("R source files: %d\n", length(r_files)))

# Count test files
test_files <- list.files("tests/testthat", pattern = "^test.*\\.R$", full.names = FALSE)
cat(sprintf("Test files: %d\n\n", length(test_files)))

# List test files
cat("Test files:\n")
for (tf in sort(test_files)) {
  # Count test_that calls in each file
  full_path <- file.path("tests/testthat", tf)
  lines <- readLines(full_path, warn = FALSE)
  n_tests <- length(grep("test_that\\(", lines))
  cat(sprintf("  %-40s  %3d tests\n", tf, n_tests))
}

total_tests <- sum(sapply(test_files, function(tf) {
  full_path <- file.path("tests/testthat", tf)
  lines <- readLines(full_path, warn = FALSE)
  length(grep("test_that\\(", lines))
}))

cat(sprintf("\nTotal test_that() calls: %d\n", total_tests))

# Estimate which R files have tests
cat("\n=== Test Coverage Estimate ===\n")
tested_modules <- gsub("^test_(.+)\\.R$", "\\1", test_files)
tested_modules <- gsub("_", "", tested_modules)  # Remove underscores for matching

r_modules <- gsub("\\.R$", "", r_files)
r_modules <- gsub("_", "", r_modules)

has_test <- tolower(r_modules) %in% tolower(tested_modules)
coverage_pct <- mean(has_test) * 100

cat(sprintf("\nFiles with tests: %d/%d (%.1f%%)\n", sum(has_test), length(r_files), coverage_pct))

cat("\nFiles WITHOUT dedicated test files:\n")
untested <- r_files[!has_test]
for (uf in sort(untested)) {
  cat(sprintf("  - %s\n", uf))
}

cat("\n=== Summary ===\n")
cat(sprintf("Estimated coverage: ~%.0f%%\n", coverage_pct))
cat("(Run 'make coverage' for detailed line-by-line coverage)\n\n")
