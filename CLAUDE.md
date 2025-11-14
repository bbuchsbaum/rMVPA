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
  - **Continuous Decoding**: hrfdecoder model for TR-level continuous-time decoding (see section below)

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
- Commonly used caret models (rf, spls, etc.) have been extracted and adapted in `R/caret_models.R`
- Model aliases exist for backward compatibility: `sda` → `sda_notune`, `glmnet` → `glmnet_opt`

## Continuous Decoding Models

### Overview

rMVPA supports **continuous-time decoding** through the `hrfdecoder` integration. Unlike traditional MVPA that operates on trial-level beta estimates, continuous decoding works directly on TR-level fMRI data.

**Key Files:**
- `R/hrfdecoder_design.R`: Specialized design object for continuous-time data
- `R/hrfdecoder_model.R`: Model adapter implementing S3 methods
- `vignettes/Continuous_Decoding.Rmd`: Comprehensive tutorial

### Architecture Pattern

The hrfdecoder integration demonstrates the **recommended pattern** for extending rMVPA with continuous decoders:

1. **Specialized Design** (`hrfdecoder_design`):
   - Extends `mvpa_design` via subclassing: `c("hrfdecoder_design", "mvpa_design", "list")`
   - Stores temporal metadata: `$event_model` (from fmridesign), `$events` (trial table)
   - Uses dummy `y_train` (TR indices) for CV fold construction only
   - Actual decoding targets come from event metadata

2. **Model Spec** (`hrfdecoder_model`):
   - Created via `create_model_spec("hrfdecoder_model", ...)`
   - Stores hyperparameters: `lambda_W`, `lambda_HRF`, `lambda_smooth`, `basis`, `window`
   - Delegates to external `hrfdecoder` package for solver

3. **S3 Method Implementation**:
   - `y_train.hrfdecoder_model()`: Returns `seq_len(nobs(dataset))` for CV
   - `train_model.hrfdecoder_model()`: Calls `hrfdecoder::hrfdecoder_fit()` on TR-level ROI data
   - `format_result.hrfdecoder_model()`: Predicts TR-level, aggregates to event-level
   - `merge_results.hrfdecoder_model()`: Builds `classification_result` from folds

### Key Differences from Standard MVPA

| Aspect | Standard MVPA | hrfdecoder |
|--------|---------------|------------|
| Data granularity | Trial-level (one obs per event) | TR-level (one obs per TR) |
| `y_train` | Event labels | Dummy TR sequence |
| `block_var` | Per-event run IDs | Per-TR run IDs |
| CV splits | Event-level | TR-level (blocked by run) |
| Training | Fits on trial betas | Fits on continuous TRs |
| Prediction | Direct classification | TR prediction → event aggregation |

### Cross-Validation

- **Fold construction**: Uses dummy `y_train` (1:T) with `block_var` (per-TR run IDs)
- **Training**: Entire runs held out (all TRs from that run)
- **Testing**: Predicts on held-out TRs, then aggregates to events via `window` parameter
- **Metrics**: Computed on event-level after aggregation (not TR-level)

### Data Flow

```
fMRI data (T × V) → mvpa_dataset
     ↓
event_model + events + block_var → hrfdecoder_design
     ↓
hrfdecoder_model(dataset, design) → model spec
     ↓
run_searchlight / run_regional → ROI extraction → train_model
     ↓                                                ↓
     |                                   hrfdecoder::hrfdecoder_fit (ALS solver)
     ↓                                                ↓
format_result → predict TR-level → aggregate_events → classification result
     ↓
merge_results → combine folds → performance metrics
```

### Dependencies

- **fmridesign**: Event model construction, basis matrices
- **fmrihrf**: HRF basis functions (spmg1, spmg3, FIR, etc.)
- **hrfdecoder**: Core continuous-time solver (C++/Armadillo)
- Runtime check via `requireNamespace("hrfdecoder", quietly = TRUE)`

### Extending for Other Continuous Decoders

To add a new continuous decoder:

1. Create `my_continuous_design()` following the hrfdecoder pattern
2. Implement S3 methods: `y_train.*`, `train_model.*`, `format_result.*`, `merge_results.*`
3. Store temporal metadata in design subclass fields
4. Use dummy `y_train` for CV, real targets from metadata
5. Aggregate continuous predictions to desired output granularity

**No core framework changes required** - the S3 system is fully extensible.

### Common Pitfalls

- **block_var length**: Must match total TRs, not number of events
- **Event timing**: Validate events don't exceed total acquisition time
- **Aggregation windows**: Events near run boundaries may have incomplete windows
- **y_train semantics**: It's a dummy variable, not the decoding target
- **Performance metrics**: Computed on aggregated events, not individual TRs

### Validation Features

The design constructor includes automatic validation:
- Checks for required columns in events table (`onset`, `condition`)
- Validates block_var length against sampling_frame
- Warns about events outside acquisition time
- Warns about events with zero probabilities after aggregation

### Reference Implementation

hrfdecoder is the **canonical example** of continuous decoding in rMVPA. Review:
- `R/hrfdecoder_design.R` for design pattern
- `R/hrfdecoder_model.R` for S3 method implementation
- `vignettes/Continuous_Decoding.Rmd` for complete workflow
- Tests in `tests/testthat/test_hrfdecoder.R` (when added)