# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

rMVPA is an R package for Multivoxel Pattern Analysis (MVPA) of neuroimaging data. It provides infrastructure for machine learning analyses on neuroimaging datasets, supporting both programmatic R usage and command-line interfaces.

## Common Development Tasks

### Building and Testing
```bash
# Install package locally
R CMD INSTALL .

# Run R CMD check
R CMD check .

# Run specific tests
Rscript -e "testthat::test_file('tests/testthat/test_mvpa_searchlight.R')"

# Run all tests
Rscript -e "testthat::test_local()"

# Build documentation
Rscript -e "devtools::document()"

# Build vignettes
Rscript -e "devtools::build_vignettes()"

# Install with dependencies
Rscript -e "devtools::install(dependencies = TRUE)"
```

### Linting and Code Quality
```bash
# Run lintr (if available)
Rscript -e "lintr::lint_package()"

# Check code coverage
Rscript -e "covr::package_coverage()"
```

## High-Level Architecture

### Core Design Patterns

1. **S3 Object System**: The package uses S3 classes extensively for model types (`mvpa_model`, `rsa_model`, `contrast_rsa_model`, etc.) with generic functions (`train_model`, `predict_model`, `performance`).

2. **Model Registry**: All machine learning models are stored in the `MVPAModels` environment (R/classifiers.R). Models follow a standard interface:
   - `type`: "Classification" or "Regression"
   - `library`: Required packages
   - `parameters`: Tunable hyperparameters
   - `grid`: Function to generate parameter grids
   - `fit`: Model training function
   - `predict`: Prediction function
   - `prob`: Probability estimation (for classification)

3. **Analysis Workflows**: Two main analysis approaches:
   - **Regional**: Analyzes specific brain regions (ROIs)
   - **Searchlight**: Sliding sphere analysis across the brain

4. **Cross-validation Infrastructure**: Custom cross-validation system in R/crossval.R supports various schemes (k-fold, blocked, bootstrap).

### Key Components

- **Data Structures**:
  - `mvpa_dataset`: Core data container with neuroimaging data, labels, and metadata
  - `mvpa_design`: Experimental design specification
  - Model-specific designs: `rsa_design`, `manova_design`, `feature_rsa_design`

- **Model Types**:
  - Standard MVPA: Classification/regression models
  - RSA (Representational Similarity Analysis): Multiple variants
  - MANOVA: Multivariate ANOVA
  - MS-ReVE: Multi-Scale Representational Variance Explained (in `contrast_rsa_model.R`)

- **Feature Selection**: Modular system supporting F-test, categorical scores, and custom methods

- **Performance Metrics**: Extensible performance evaluation supporting accuracy, AUC, RMSE, R², and custom metrics

### Important Ongoing Work

**Caret Removal**: The package is undergoing a major refactoring to remove the `caret` dependency (see `caret_removal_plan.md`). Key changes:
- Replacing `caret::train` with custom tuning loops using `rsample` and `yardstick`
- Model loading now exclusively uses the internal `MVPAModels` registry
- Performance metrics transitioning from `caret` to `yardstick` functions

Current branch `remaining-codex-integration` contains work-in-progress changes.

### Command Line Scripts

The `scripts/` directory contains CLI tools:
- `MVPA_Searchlight.R`: Searchlight analysis
- `MVPA_Regional.R`: Regional/ROI analysis
- `MVPA_Cluster.R`: Cluster-based analysis
- `MVPA_Predict.R`: Prediction on new data

### Testing Infrastructure

- Comprehensive test suite in `tests/testthat/`
- Tests organized by functionality (searchlight, regional, RSA, etc.)
- Mock CV specifications helper in `helper-mock_cv_spec.R`
- GitHub Actions CI runs R CMD check on multiple platforms

### Key Dependencies

- **Neuroimaging**: `neuroim2`, `neurosurf` (custom packages by the same author)
- **Machine Learning**: `rsample` (resampling), `yardstick` (metrics) from tidymodels ecosystem
- **Data Manipulation**: `dplyr`, `purrr`, `tibble`
- **Parallel Processing**: `future`, `future.apply`, `furrr`
- **Statistical Models**: `glmnet`, `sda`, `randomForest`, `e1071`

### Model Registry (MVPAModels)

The package uses a centralized model registry called `MVPAModels` (an environment in `R/classifiers.R`). All built-in models are registered here with a consistent structure:

- **Model specification fields**:
  - `type`: "Classification" or "Regression"
  - `library`: Required R packages for the model
  - `label`: Display name
  - `parameters`: data.frame of tunable parameters
  - `grid`: Function to generate parameter tuning grid
  - `fit`: Function to train the model
  - `predict`: Function for predictions
  - `prob`: Function for class probabilities (classification only)

- **Loading models**: Use `load_model(name)` to retrieve a model specification
- **Custom models**: Use `register_mvpa_model(name, spec)` to add new models

### Development Notes

- The package heavily uses functional programming patterns with `purrr`
- Parallel processing is handled through the `future` framework
- Logging via `futile.logger` for debugging
- MS-ReVE functionality is in `R/contrast_rsa_model.R` (not `contrast_rsa.R`)
- The package no longer depends on `caret` - all model management is handled internally