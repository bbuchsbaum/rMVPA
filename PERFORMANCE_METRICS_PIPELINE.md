# rMVPA Performance Metrics Pipeline: Complete Data Flow Analysis

## Executive Summary

This document traces the complete pipeline from model predictions to final performance metrics in rMVPA's searchlight and regional analyses. It identifies the required data structures at each stage and compares hrfdecoder's current implementation against the standard contracts.

---

## 1. Pipeline Entry Point: `mvpa_iterate()`

**Location**: `R/mvpa_iterate.R` (line 366)

### Flow:
1. Iterates over ROI batches
2. Calls `extract_roi()` to get ROI data
3. Calls `process_roi()` dispatcher for each ROI
4. Returns tibble with columns: `result`, `indices`, `performance`, `id`, `error`, `error_message`

**Key insight**: `process_roi()` is a generic dispatcher that routes to `internal_crossval()` or `external_crossval()` based on whether the model has a test set.

---

## 2. Cross-Validation: `internal_crossval()`

**Location**: `R/mvpa_iterate.R` (line 224)

### Flow:
```
1. Generate CV samples: crossval_samples(mspec$crossval, data, y_train(mspec))
   ├─ For hrfdecoder: y_train = seq_len(nobs(dataset))  [DUMMY TR INDICES]
   └─ For mvpa_model: y_train = actual labels

2. For each fold (via pmap):
   ├─ Call train_model(mspec, train, ytrain, sl_info, cv_spec, indices)
   └─ Call format_result(mspec, result, error_message, context)

3. Bind all fold results → result_set (tibble)

4. Call merge_results(mspec, result_set, indices, id)
   └─ Returns final tibble with performance
```

### Critical Context Object:
```r
context = list(
  roi = roi,
  ytrain = ytrain,  # Fold-specific training labels
  ytest = ytest,    # Fold-specific test labels
  train = train,    # Training data (resample object)
  test = test,      # Test data (resample object)
  .id = .id         # Fold identifier
)
```

---

## 3. Per-Fold Processing: `train_model()` → `format_result()`

### 3.1 Standard MVPA Model (`mvpa_model`)

**train_model** (`R/mvpa_model.R`, called from various model implementations):
- Receives: `(obj, train_dat, ytrain, sl_info, cv_spec, indices)`
- Returns: Fitted model object

**format_result.mvpa_model** (`R/mvpa_model.R:81`):
```r
format_result.mvpa_model <- function(obj, result, error_message=NULL, context, ...) {
  if (!is.null(error_message)) {
    return(tibble(class=list(NULL), probs=list(NULL), y_true=list(context$ytest),
                  fit=list(NULL), error=TRUE, error_message=error_message))
  }

  # Predict on test set
  pred <- predict(result, tibble::as_tibble(context$test), NULL)

  # Build output
  plist <- lapply(pred, list)
  plist$y_true <- list(context$ytest)
  plist$test_ind <- list(as.integer(context$test))
  plist$fit <- if (obj$return_fit) list(result) else list(NULL)
  plist$error <- FALSE
  plist$error_message <- "~"

  tibble::as_tibble(plist)
}
```

**Required Output Structure**:
```r
tibble(
  class = list(predicted_classes),    # Factor vector of predictions
  probs = list(probability_matrix),   # Matrix: rows=obs, cols=classes
  y_true = list(ytest),               # True labels from context$ytest
  test_ind = list(test_indices),      # Integer indices of test observations
  fit = list(model_fit_or_NULL),      # Optional fitted model
  error = FALSE,
  error_message = "~"
)
```

### 3.2 hrfdecoder Model

**train_model.hrfdecoder_model** (`R/hrfdecoder_model.R:193`):
- **IGNORES `y` parameter** (it's a dummy TR sequence)
- Uses `obj$design$event_model` and `obj$design$events` for targets
- Calls `hrfdecoder::fit_hrfdecoder()` on TR-level data
- Returns: `hrfdecoder_fit_wrap` object

**format_result.hrfdecoder_model** (`R/hrfdecoder_model.R:239`):
```r
format_result.hrfdecoder_model <- function(obj, result, error_message=NULL, context, ...) {
  # 1. Predict TR-level soft labels on test TRs
  Ptest <- hrfdecoder::predict_hrfdecoder(object=result$fit, Y_test=Xtest, mode="tr")

  # 2. Aggregate TR predictions to event-level probabilities
  agg <- hrfdecoder::aggregate_events(
    P = Ptest,
    events = obj$design$events,
    TR = result$fit$settings$TR,
    window = obj$window,
    hrf = result$fit$hrf,
    normalize = TRUE
  )

  # 3. Extract event-level outputs
  probs <- as.matrix(agg$probs)        # Events × Classes
  observed <- agg$y_true               # Event labels

  # 4. Return in standard format
  plist <- list(
    class = list(pred_class),
    probs = list(probs),
    y_true = list(observed),           # EVENT-level labels (from aggregation)
    test_ind = list(as.integer(context$test)),  # TR indices (ignored downstream)
    fit = list(...),
    error = FALSE,
    error_message = "~"
  )
  tibble::as_tibble(plist)
}
```

**Key Differences**:
- `y_true` comes from `aggregate_events()`, NOT from `context$ytest`
- `context$ytest` contains dummy TR indices (1:T), which are IGNORED
- `test_ind` points to TRs, but downstream only uses `y_true` and predictions

---

## 4. Fold Merging: `merge_results()`

### 4.1 Standard MVPA Model

**merge_results.mvpa_model** (`R/mvpa_model.R:50`):
```r
merge_results.mvpa_model <- function(obj, result_set, indices, id, ...) {
  # Check for errors
  if (any(result_set$error)) { return error tibble }

  # Wrap results across folds
  cres <- if (obj$return_fit) {
    predictor <- weighted_model(result_set$fit)
    wrap_result(result_set, obj$design, predictor)
  } else {
    wrap_result(result_set, obj$design)
  }

  # Compute performance
  if (obj$compute_performance) {
    tibble(result=list(cres), indices=list(indices),
           performance=list(compute_performance(obj, cres)),
           id=id, error=FALSE, error_message="~")
  } else {
    tibble(result=list(cres), indices=list(indices),
           performance=list(NULL), id=id, error=FALSE, error_message="~")
  }
}
```

**wrap_result()** (`R/mvpa_model.R:3`):
- Aggregates probabilities across folds for each test observation
- Calls `classification_result()` or `regression_result()` constructor
- Returns: Classification/regression result object

### 4.2 hrfdecoder Model

**merge_results.hrfdecoder_model** (`R/hrfdecoder_model.R:346`):
```r
merge_results.hrfdecoder_model <- function(obj, result_set, indices, id, ...) {
  if (any(result_set$error)) { return error tibble }

  # Concatenate per-fold event-level outputs
  probs_list <- lapply(result_set$probs, as.matrix)
  probs <- do.call(rbind, probs_list)      # Stack all event predictions
  y_true <- do.call(c, result_set$y_true)  # Concatenate all event labels
  predicted <- do.call(c, result_set$class)

  # Build standard classification_result
  cres <- classification_result(y_true, predicted, probs,
                                testind=NULL, test_design=NULL, predictor=NULL)

  # Compute performance
  perf <- if (obj$compute_performance) compute_performance(obj, cres) else NULL

  tibble(result=list(cres), indices=list(indices),
         performance=list(perf), id=id, error=FALSE, error_message="~")
}
```

**Key Differences**:
- Does NOT use `wrap_result()` (which expects trial-aligned test indices)
- Directly concatenates event-level predictions from all folds
- Each fold contributes different events (blocked by run)
- No probability averaging needed (each event appears once)

---

## 5. Result Constructor: `classification_result()`

**Location**: `R/mvpa_result.R:35`

```r
classification_result <- function(observed, predicted, probs,
                                  testind=NULL, test_design=NULL, predictor=NULL) {
  if (is.numeric(observed)) {
    regression_result(observed, predicted, testind, test_design, predictor)
  } else if (length(levels(as.factor(observed))) == 2) {
    binary_classification_result(as.factor(observed), predicted, probs, ...)
  } else if (length(levels(as.factor(observed))) > 2) {
    multiway_classification_result(as.factor(observed), predicted, probs, ...)
  }
}
```

### Required Fields:
```r
structure(list(
  observed = factor(...),      # True class labels
  predicted = factor(...),     # Predicted class labels
  probs = matrix(...),         # Probability matrix (obs × classes)
  testind = integer_or_NULL,   # Optional test indices
  test_design = df_or_NULL,    # Optional test design
  predictor = obj_or_NULL      # Optional predictor object
), class = c("multiway_classification_result", "classification_result", "list"))
```

**Critical**: Only `observed`, `predicted`, and `probs` are used for performance computation. The `testind` field is optional and only used for subsetting in some contexts.

---

## 6. Performance Computation: `compute_performance()`

**Generic**: `compute_performance(obj, result)`

### 6.1 Standard Models

**compute_performance.mvpa_model** (`R/mvpa_model.R:133`):
```r
compute_performance.mvpa_model <- function(obj, result) {
  obj$performance(result)
}
```

Where `obj$performance` is set during model construction to one of:
- `get_binary_perf(split_list)` → calls `performance.binary_classification_result()`
- `get_multiclass_perf(split_list, class_metrics)` → calls `performance.multiway_classification_result()`
- `get_regression_perf(split_list)` → calls `performance.regression_result()`

### 6.2 hrfdecoder Model

**compute_performance.hrfdecoder_model** (`R/hrfdecoder_model.R:377`):
```r
compute_performance.hrfdecoder_model <- function(obj, result) {
  obj$performance(result)
}
```

**Default**: `get_multiclass_perf(design$split_groups, class_metrics=TRUE)` set at line 143.

---

## 7. Performance Metrics: `performance()` methods

**Location**: `R/performance.R`

### Binary Classification

**performance.binary_classification_result** (`R/performance.R:172`):
```r
performance.binary_classification_result <- function(x, split_list=NULL, ...) {
  if (is.null(split_list)) {
    ret <- binary_perf(x$observed, x$predicted, x$probs)
  } else {
    # Compute overall + per-split metrics
  }
}
```

**binary_perf()** (`R/performance.R:234`):
```r
binary_perf <- function(observed, predicted, probs) {
  res_acc <- yardstick::accuracy_vec(truth=observed, estimate=predicted)
  prob_positive <- probs[, levels(observed)[2]]
  res_auc <- yardstick::roc_auc_vec(truth=observed, estimate=prob_positive,
                                    event_level="second")
  c(Accuracy=res_acc, AUC=res_auc - 0.5)  # Note: subtracts 0.5 from AUC
}
```

### Multiway Classification

**performance.multiway_classification_result** (`R/performance.R:196`):
```r
performance.multiway_classification_result <- function(x, split_list=NULL,
                                                       class_metrics=FALSE, ...) {
  multiclass_perf(x$observed, x$predicted, x$probs, class_metrics)
}
```

**multiclass_perf()** (`R/performance.R:259`):
```r
multiclass_perf <- function(observed, predicted, probs, class_metrics=FALSE) {
  acc <- yardstick::accuracy_vec(truth=observed, estimate=predicted)

  # Per-class AUC (one-vs-rest)
  aucres <- sapply(seq_along(levels(observed)), function(i) {
    lev <- levels(observed)[i]
    pos <- observed == lev
    pclass <- probs[,i]
    pother <- rowMeans(probs[,-i, drop=FALSE])
    score <- pclass - pother
    binary_truth <- factor(ifelse(pos, "positive", "negative"),
                          levels=c("negative", "positive"))
    yardstick::roc_auc_vec(truth=binary_truth, estimate=score,
                          event_level="second") - 0.5
  })

  names(aucres) <- paste0("AUC_", colnames(probs))
  metrics <- c(Accuracy=acc, AUC=mean(aucres, na.rm=TRUE))

  if (class_metrics) c(metrics, aucres) else metrics
}
```

### Regression

**performance.regression_result** (`R/performance.R:39`):
```r
performance.regression_result <- function(x, split_list=NULL, ...) {
  obs <- x$observed
  pred <- x$predicted

  res_rsq <- yardstick::rsq_vec(truth=obs, estimate=pred)
  res_rmse <- yardstick::rmse_vec(truth=obs, estimate=pred)
  res_spearcor <- cor(obs, pred, method="spearman", use="pairwise.complete.obs")

  c(R2=res_rsq, RMSE=res_rmse, spearcor=res_spearcor)
}
```

---

## 8. Complete Data Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────────┐
│                         mvpa_iterate()                                  │
│  Entry point for searchlight/regional analysis                         │
└──────────────────────┬──────────────────────────────────────────────────┘
                       │
                       ▼
         ┌─────────────────────────┐
         │  process_roi.default()  │
         │  Dispatcher              │
         └────────┬────────────────┘
                  │
                  ├─ has_test_set? → external_crossval()
                  └─ else → internal_crossval()
                              │
                              ▼
        ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
        ┃       internal_crossval()              ┃
        ┃  1. Generate CV samples                ┃
        ┃  2. For each fold:                     ┃
        ┗━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┛
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼                         ▼
   train_model()           context = {ytrain, ytest, train, test}
        │
        └─────────────┬───────────┘
                      │
                      ▼
              format_result()
        ┌─────────────┴─────────────┐
        │                           │
        │  Standard MVPA            │  hrfdecoder
        │  ─────────────            │  ──────────
        │  1. Predict on test       │  1. Predict TR-level (test TRs)
        │  2. Use context$ytest     │  2. Aggregate to events
        │     as y_true             │  3. Use agg$y_true (NOT context$ytest)
        │                           │  4. context$ytest = dummy TRs (IGNORED)
        └───────────┬───────────────┘
                    │
        ┌───────────▼────────────┐
        │  Per-fold output:      │
        │  tibble(               │
        │    class = list(...),  │
        │    probs = list(...),  │
        │    y_true = list(...), │
        │    test_ind = list(...),
        │    fit = list(...),    │
        │    error = FALSE       │
        │  )                     │
        └────────┬───────────────┘
                 │
                 ▼
        bind_rows() → result_set
                 │
                 ▼
        ┏━━━━━━━━━━━━━━━━━━━━━┓
        ┃  merge_results()    ┃
        ┗━━━━━━━━┳━━━━━━━━━━━┛
                 │
    ┌────────────┴────────────┐
    │                         │
    │  Standard MVPA          │  hrfdecoder
    │  ─────────────          │  ──────────
    │  wrap_result()          │  Direct concatenation:
    │  ├─ Aggregate probs     │  ├─ rbind(probs)
    │  │  across folds        │  ├─ c(y_true)
    │  ├─ Use design$y_test   │  └─ c(predicted)
    │  └─ classification_     │
    │     result()            │  classification_result()
    └──────────┬──────────────┘
               │
               ▼
    ┌─────────────────────┐
    │ classification_     │
    │ result object       │
    │ ─────────────       │
    │ $observed           │ ← From y_true field
    │ $predicted          │ ← From class field
    │ $probs              │ ← Probability matrix
    │ $testind (optional) │
    └──────────┬──────────┘
               │
               ▼
    compute_performance(obj, result)
               │
               ▼
    obj$performance(result)  ← Function closure set during model construction
               │
               ├─ Binary: get_binary_perf() → binary_perf()
               ├─ Multiclass: get_multiclass_perf() → multiclass_perf()
               └─ Regression: get_regression_perf()
               │
               ▼
    ┌─────────────────────────┐
    │ Performance Metrics     │
    │ ─────────────────────   │
    │ Binary:                 │
    │   - Accuracy            │
    │   - AUC                 │
    │                         │
    │ Multiclass:             │
    │   - Accuracy            │
    │   - AUC (mean)          │
    │   - AUC per class       │
    │                         │
    │ Regression:             │
    │   - R²                  │
    │   - RMSE                │
    │   - Spearman cor        │
    └─────────────────────────┘
               │
               ▼
    Return tibble(result, indices, performance, id, error, error_message)
```

---

## 9. Contract Requirements Summary

### `format_result()` Must Return:

```r
tibble(
  class = list(factor_vector),        # Predicted classes
  probs = list(matrix),               # Probability matrix (obs × classes)
  y_true = list(factor_or_numeric),   # TRUE LABELS FOR OBSERVATIONS
  test_ind = list(integer_vector),    # Test observation indices (optional)
  fit = list(model_or_NULL),          # Optional fitted model
  error = logical(1),                 # Error flag
  error_message = character(1)        # Error message or "~"
)
```

**Critical**: `y_true` MUST contain the actual target labels for the observations being predicted. This is used directly by `merge_results()` and `classification_result()`.

### `merge_results()` Must Return:

```r
tibble(
  result = list(classification_result_object),
  indices = list(roi_indices),
  performance = list(named_numeric_vector_or_NULL),
  id = scalar,
  error = logical(1),
  error_message = character(1)
)
```

### `classification_result` Must Contain:

```r
list(
  observed = factor_vector,      # Ground truth labels
  predicted = factor_vector,     # Predicted labels (same length as observed)
  probs = matrix,                # Probability matrix (nrow=length(observed))
  testind = integer_or_NULL,     # Optional
  test_design = df_or_NULL,      # Optional
  predictor = object_or_NULL     # Optional
)
```

---

## 10. hrfdecoder Compliance Analysis

### ✅ **COMPLIANT**: hrfdecoder correctly implements the pipeline

#### Evidence:

1. **format_result.hrfdecoder_model** (line 239):
   - Returns all required fields: `class`, `probs`, `y_true`, `test_ind`, `fit`, `error`, `error_message`
   - ✅ `y_true` contains **event-level labels** from `aggregate_events()`, NOT dummy TRs
   - ✅ Correctly ignores `context$ytest` (which contains dummy TR indices)

2. **merge_results.hrfdecoder_model** (line 346):
   - Concatenates event-level predictions across folds
   - Calls `classification_result()` with proper structure
   - Returns tibble with `result`, `indices`, `performance`, `id`, `error`, `error_message`
   - ✅ No dependency on trial-level alignment

3. **classification_result construction**:
   - `y_true`, `predicted`, and `probs` all have event-level granularity
   - Probabilities are normalized per-event
   - ✅ Passes directly to `multiclass_perf()` → `yardstick` functions

4. **Performance computation**:
   - Uses standard `compute_performance.hrfdecoder_model()` → `obj$performance(result)`
   - Defaults to `get_multiclass_perf()` with `class_metrics=TRUE`
   - ✅ Works identically to standard MVPA models

### Key Architectural Decision:

**The dummy `y_train` exists ONLY for CV fold construction, NOT for predictions.**

- Cross-validation splits TRs by run (via `block_var`)
- Each fold trains on N-1 runs of TR-level data
- Predictions happen at TR-level, then aggregate to events
- Aggregated event labels (from `aggregate_events()`) flow to performance metrics
- `context$ytest` (dummy TR indices) is **intentionally ignored** in `format_result`

### No Pipeline Breaks:

❌ **NO** breaks in the data flow
❌ **NO** missing transformations
❌ **NO** structural mismatches
✅ **ALL** contracts satisfied
✅ **FULL** compatibility with `yardstick` metrics

---

## 11. Comparison: Standard MVPA vs hrfdecoder

| Aspect | Standard MVPA | hrfdecoder |
|--------|---------------|------------|
| **Data Granularity** | Trial-level (1 obs per event) | TR-level (1 obs per TR) |
| **y_train Source** | `design$y_train` (event labels) | `seq_len(nobs)` (dummy TR indices) |
| **CV Splits** | Event-level (blocked by run) | TR-level (blocked by run) |
| **Training** | Fits on trial betas | Fits on continuous TRs |
| **Prediction** | Direct trial classification | TR prediction → event aggregation |
| **y_true in format_result** | `context$ytest` (trial labels) | `aggregate_events()$y_true` (event labels) |
| **context$ytest Usage** | Used as ground truth | **IGNORED** (contains dummy TRs) |
| **merge_results Strategy** | `wrap_result()` (averages probs) | Direct concatenation (no averaging) |
| **classification_result** | Trial-aligned observations | Event-aligned observations |
| **Performance Metrics** | Computed on trials | Computed on events |

### Critical Insight:

**hrfdecoder treats `context$ytest` as opaque because it contains dummy values that are meaningless for performance evaluation.** The real labels come from the event metadata stored in the design object, which is accessed via `aggregate_events()`.

This is **by design** and **fully supported** by the S3 dispatch system. The `format_result()` contract only requires that `y_true` contain meaningful labels—it doesn't mandate where those labels come from.

---

## 12. Validation Checklist

### hrfdecoder Implementation:

- [x] `y_train.hrfdecoder_model()` returns TR-count-length dummy sequence
- [x] `train_model.hrfdecoder_model()` ignores `y` parameter
- [x] `format_result.hrfdecoder_model()` aggregates TR predictions to events
- [x] Event labels extracted from `aggregate_events()`, not `context$ytest`
- [x] Returns tibble with all required fields (`class`, `probs`, `y_true`, ...)
- [x] `merge_results.hrfdecoder_model()` concatenates event-level outputs
- [x] Calls `classification_result()` with event-level data
- [x] `compute_performance.hrfdecoder_model()` uses standard performance function
- [x] Compatible with `yardstick::accuracy_vec()` and `yardstick::roc_auc_vec()`
- [x] No structural mismatches in result objects
- [x] All contracts satisfied at every pipeline stage

### Potential Edge Cases:

⚠️ **Events with zero probabilities**: Already handled with warning (line 319-323)
⚠️ **Events outside TR range**: Validated in design constructor (line 89-95)
⚠️ **Single run**: Warning issued during model construction (line 118-120)

---

## 13. Conclusion

### Summary:

The hrfdecoder integration is **fully compliant** with rMVPA's performance metrics pipeline. The key architectural pattern—using dummy `y_train` for fold construction while pulling real labels from design metadata—is elegant and requires **no modifications to the core framework**.

### Data Flow Correctness:

1. ✅ Dummy `y_train` used only for CV fold creation
2. ✅ TR-level predictions aggregate to event-level within `format_result()`
3. ✅ Event labels flow from `aggregate_events()` → `y_true` → `classification_result`
4. ✅ Standard `performance()` methods work without modification
5. ✅ `yardstick` metrics receive properly structured inputs

### Design Strengths:

- **S3 extensibility**: No changes to `internal_crossval()`, `merge_results()`, or `performance()`
- **Clean separation**: Event metadata in design, not polluting CV machinery
- **Type safety**: `classification_result()` enforces contracts at construction
- **Reusability**: Pattern generalizes to any continuous decoder

### No Breaks Identified:

After tracing the complete pipeline from `mvpa_iterate()` → `internal_crossval()` → `format_result()` → `merge_results()` → `classification_result()` → `compute_performance()` → `performance()` → `yardstick` metrics, **no data flow breaks were found**.

The hrfdecoder implementation follows the established contracts at every stage and produces performance metrics that are structurally identical to standard MVPA models.

---

## 14. Files Referenced

- `/Users/bbuchsbaum/code/rMVPA/R/mvpa_iterate.R` - Main iteration logic
- `/Users/bbuchsbaum/code/rMVPA/R/crossval.R` - CV infrastructure
- `/Users/bbuchsbaum/code/rMVPA/R/mvpa_model.R` - Standard model reference
- `/Users/bbuchsbaum/code/rMVPA/R/mvpa_result.R` - Result constructors
- `/Users/bbuchsbaum/code/rMVPA/R/performance.R` - Metric computation
- `/Users/bbuchsbaum/code/rMVPA/R/hrfdecoder_model.R` - hrfdecoder adapter
- `/Users/bbuchsbaum/code/rMVPA/R/hrfdecoder_design.R` - Design specialization
- `/Users/bbuchsbaum/code/rMVPA/R/allgeneric.R` - Generic definitions

---

**Generated**: 2025-11-10
**Package**: rMVPA
**Analysis**: Complete performance metrics pipeline trace with hrfdecoder compliance verification
