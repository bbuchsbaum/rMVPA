# Makefile for rMVPA R package development
#
# Common targets:
#   make check      - Run R CMD check
#   make test       - Run all tests
#   make coverage   - Generate test coverage report
#   make document   - Update documentation
#   make install    - Install package locally
#   make clean      - Clean build artifacts

.PHONY: all check test test-perf coverage document install build clean help

# Default target
all: document check

# Run R CMD check
check:
	@echo "Running R CMD check..."
	Rscript -e "devtools::check()"

# Run tests only (faster than full check)
test:
	@echo "Running tests..."
	RGL_USE_NULL=TRUE Rscript -e "devtools::test()"

# Run performance guardrail tests
test-perf:
	@echo "Running performance guardrail tests..."
	RGL_USE_NULL=TRUE RMVPA_RUN_PERF_TESTS=true Rscript scripts/run_perf_guardrails.R

# Generate test coverage report
coverage:
	@echo "Generating test coverage report..."
	@Rscript -e "if (!requireNamespace('covr', quietly = TRUE)) install.packages('covr')"
	@Rscript -e "devtools::load_all(); cov <- covr::package_coverage(type = 'tests', quiet = FALSE); print(cov); covr::report(cov, file = 'coverage.html'); cat('\nCoverage report written to coverage.html\n')"

# Quick coverage summary (no HTML report)
coverage-summary:
	@echo "Computing test coverage (this may take several minutes)..."
	@Rscript scripts/quick_test_summary.R

# Full coverage with covr (slow, generates HTML report)
coverage-full:
	@echo "Running full coverage analysis (this may take 10+ minutes)..."
	@Rscript -e "if (!requireNamespace('covr', quietly = TRUE)) install.packages('covr')"
	@Rscript -e "devtools::load_all(quiet = TRUE); cov <- covr::package_coverage(type = 'tests', quiet = FALSE); print(cov); cat(sprintf('\n\nOverall coverage: %.1f%%\n', covr::percent_coverage(cov)))"

# Update documentation (roxygen2)
document:
	@echo "Updating documentation..."
	Rscript -e "devtools::document()"

# Install package locally
install:
	@echo "Installing package..."
	Rscript -e "devtools::install()"

# Build package tarball
build:
	@echo "Building package..."
	Rscript -e "devtools::build()"

# Build vignettes
vignettes:
	@echo "Building vignettes..."
	Rscript -e "devtools::build_vignettes()"

# Build pkgdown site
site:
	@echo "Building pkgdown site..."
	Rscript -e "pkgdown::build_site()"

# Run lintr
lint:
	@echo "Running lintr..."
	@Rscript -e "if (!requireNamespace('lintr', quietly = TRUE)) install.packages('lintr')"
	@Rscript -e "lintr::lint_package()"

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf rMVPA.Rcheck/
	@rm -f rMVPA_*.tar.gz
	@rm -f coverage.html
	@rm -rf docs/
	@find . -name "*.o" -delete
	@find . -name "*.so" -delete
	@find . -name "*.dll" -delete
	@echo "Clean complete"

# Help target
help:
	@echo "rMVPA Makefile targets:"
	@echo ""
	@echo "  make check            - Run R CMD check"
	@echo "  make test             - Run all tests"
	@echo "  make test-perf        - Run performance guardrail tests"
	@echo "  make coverage         - Generate HTML coverage report (slow)"
	@echo "  make coverage-summary - Quick test file count summary (fast)"
	@echo "  make coverage-full    - Full line-by-line coverage analysis (very slow)"
	@echo "  make document         - Update roxygen2 documentation"
	@echo "  make install          - Install package locally"
	@echo "  make build            - Build package tarball"
	@echo "  make vignettes        - Build vignettes"
	@echo "  make site             - Build pkgdown site"
	@echo "  make lint             - Run lintr checks"
	@echo "  make clean            - Clean build artifacts"
	@echo "  make help             - Show this help message"
	@echo ""
	@echo "Quick workflow:"
	@echo "  make document && make test    - Update docs and test"
	@echo "  make coverage-summary         - Quick coverage estimate"
	@echo "  make check                    - Full CRAN check"
