# Continuous Decoding Integration - Summary

## Overview

This document summarizes the improvements made to the hrfdecoder continuous decoding integration in rMVPA.

## Completed Work

### 1. Code Consolidation ✅
- **Issue**: Duplicate `hrfdecoder_model.R` implementations existed in both rMVPA and hrfdecoder repositories
- **Resolution**: Confirmed rMVPA version as canonical (better error handling, uses `classification_result`)
- **Location**: `/Users/bbuchsbaum/code/rMVPA/R/hrfdecoder_model.R`

### 2. Enhanced Validation in `hrfdecoder_design()` ✅
**File**: `R/hrfdecoder_design.R`

Added validation checks:
- Verifies events table contains required columns (`onset`, `condition`)
- Validates `block_var` length matches `sampling_frame` total TRs
- Warns if events have onsets beyond total acquisition time
- Improved roxygen2 documentation with full example

### 3. Improved Error Handling in `hrfdecoder_model()` ✅
**File**: `R/hrfdecoder_model.R`

Enhancements:
- Added comprehensive parameter documentation
- Full workflow example in roxygen2 `@examples`
- Better explanation of continuous-time decoding approach
- Warning for events with zero probabilities (empty aggregation windows)

### 4. New `print.hrfdecoder_design()` Method ✅
**File**: `R/hrfdecoder_design.R`

Pretty-prints design objects showing:
- Number of TRs
- Number of runs/blocks
- Event count and conditions
- TR duration and block lengths
- Split group information

### 5. Comprehensive Vignette ✅
**File**: `vignettes/Continuous_Decoding.Rmd`

Complete tutorial covering:
- Introduction to continuous vs trial-level MVPA
- Step-by-step workflow (event table → searchlight)
- Cross-validation strategy explanation
- Advanced topics (custom metrics, baselines, edge cases)
- Guide for extending rMVPA with other continuous decoders
- Troubleshooting section

### 6. Integration Tests ✅
**File**: `tests/testthat/test_hrfdecoder.R`

Test coverage:
- Event table validation (missing columns)
- block_var length validation
- Events outside acquisition time warning
- Proper class hierarchy (`hrfdecoder_design` → `mvpa_design`)
- Print method functionality
- Model spec creation
- `y_train()` returns TR sequence
- Design type checking in model constructor

All tests skip gracefully if hrfdecoder not installed.

### 7. Updated CLAUDE.md Documentation ✅
**File**: `CLAUDE.md`

Added comprehensive "Continuous Decoding Models" section covering:
- Architecture pattern (design + model + S3 methods)
- Key differences from standard MVPA
- Cross-validation strategy
- Data flow diagram
- Dependencies
- Extension guide for new continuous decoders
- Common pitfalls and validation features

## Key Design Principles

### 1. No Core Framework Changes
The existing rMVPA S3 infrastructure already supports continuous decoding:
- `mvpa_dataset` naturally handles TR×voxel matrices
- `mvpa_design` subclassing pattern accommodates temporal metadata
- Generic methods (`train_model.*`, `format_result.*`, `merge_results.*`) enable full customization
- CV machinery is observation-agnostic

### 2. Clean Separation of Concerns
- **hrfdecoder package**: Pure solver (C++/Armadillo optimization)
- **rMVPA integration**: Model adapter (design + S3 methods)
- **fmridesign/fmrihrf**: Temporal modeling infrastructure

### 3. Pattern for Future Extensions
hrfdecoder serves as the **canonical reference implementation** for:
- Continuous-time decoders
- State-space models
- Temporal CNN/RNN decoders
- Any model operating on continuous observations

## Files Modified

```
R/hrfdecoder_design.R              # Enhanced validation + print method + examples
R/hrfdecoder_model.R               # Better docs + error handling + examples
vignettes/Continuous_Decoding.Rmd  # Comprehensive tutorial (NEW)
tests/testthat/test_hrfdecoder.R   # Integration tests (NEW)
CLAUDE.md                          # Continuous decoding section
CONTINUOUS_DECODING_SUMMARY.md     # This file (NEW)
```

## Usage Example

```r
library(rMVPA)
library(hrfdecoder)
library(fmridesign)
library(fmrihrf)

# 1. Event table
events_df <- data.frame(
  onset = seq(10, 290, by = 20),
  condition = rep(c("A", "B", "C"), length.out = 15),
  run = rep(1:3, each = 5)
)

# 2. Sampling frame
sframe <- sampling_frame(blocklens = c(100, 100, 100), TR = 2)

# 3. Event model
evmod <- event_model(
  onset ~ hrf(condition, basis = "spmg3"),
  data = events_df,
  block = ~run,
  sampling_frame = sframe
)

# 4. rMVPA structures
dset <- mvpa_dataset(train_data = fmri_data, mask = mask)

design <- hrfdecoder_design(
  event_model = evmod,
  events = events_df,
  block_var = rep(1:3, each = 100)
)

# 5. Model specification
mspec <- hrfdecoder_model(
  dataset = dset,
  design = design,
  lambda_W = 10,
  lambda_HRF = 1,
  lambda_smooth = 5,
  basis = fmrihrf::spmg1(),
  window = c(4, 8)
)

# 6. Analysis
results <- run_searchlight(mspec, radius = 8, method = "randomized")
```

## Next Steps (Future Work)

Optional enhancements not included in this minimal approach:

1. **Hyperparameter Tuning**: Add nested CV support for `lambda_W`, `lambda_HRF`, `lambda_smooth`
2. **TR-Level Outputs**: Option to return TR-level soft labels (not just event-aggregated)
3. **Continuous Regression**: Extend to continuous response variables
4. **Multi-Subject Support**: Subject-level blocking and cross-subject aggregation
5. **Performance Optimizations**: Parallel class updates, prewhitening, caching

## Testing

Run tests:
```bash
Rscript -e "testthat::test_file('tests/testthat/test_hrfdecoder.R')"
```

Build documentation:
```bash
Rscript -e "devtools::document()"
```

Build vignettes:
```bash
Rscript -e "devtools::build_vignettes()"
```

## References

- **rMVPA**: `/Users/bbuchsbaum/code/rMVPA`
- **hrfdecoder**: `~/code/hrfdecoder`
- **fmridesign**: Required for event modeling
- **fmrihrf**: Required for HRF basis functions
