# Performance Analysis: rMVPA Core Engine

This document identifies performance bottlenecks and optimization opportunities in rMVPA's core inner loops.

## Executive Summary

The main performance-critical paths are:
1. **ROI extraction and filtering** (`extract_roi`, `filter_roi`) - called once per searchlight center
2. **Cross-validation sampling** (`crossval_samples`, `internal_crossval`) - called once per ROI
3. **Results combination** (`combine_randomized`, `combine_standard`) - called once per searchlight run

For a typical searchlight analysis with 50,000 centers and 5 CV folds, these functions are called millions of times.

---

## Critical Performance Issues

### 1. `rowwise()` + `mutate()` Pattern in ROI Extraction

**Location**: `R/mvpa_iterate.R:440-447`

```r
sf <- sf %>%
  rowwise() %>%
  mutate(roi=list(extract_roi(sample, dset, center_global_id = rnum, min_voxels = min_voxels_required))) %>%
  select(-sample)
```

**Problem**: `dplyr::rowwise()` combined with `mutate()` has significant overhead due to:
- Creating grouped tibble structure
- Row-by-row evaluation with NSE overhead
- Method dispatch for each row

**Impact**: For 10,000 ROIs per batch, this adds ~0.5-2s overhead per batch.

**Recommended Fix**:
```r
# Use lapply or purrr::map2 directly
rois <- lapply(seq_len(nrow(sf)), function(i) {
  extract_roi(sf$sample[[i]], dset,
              center_global_id = sf$rnum[i],
              min_voxels = min_voxels_required)
})
sf$roi <- rois
sf$sample <- NULL
```

---

### 2. Per-Column `apply()` in `filter_roi.ROIVec`

**Location**: `R/resample.R:76-79`

```r
nas <- apply(trdat, 2, function(v) any(is.na(v)))
sdnonzero <- apply(trdat, 2, sd, na.rm=TRUE) > 0
```

**Problem**: `apply()` with custom functions is slow because:
- It coerces to matrix then applies function row/column-wise
- For SD calculation, it's computing full SD when we only need to know if > 0

**Impact**: Called for every searchlight center. With 50,000 centers and ~100 voxels per ROI, this is 5M+ column evaluations.

**Recommended Fix**:
```r
# For NA check - vectorized
nas <- colSums(is.na(trdat)) > 0

# For SD check - use matrixStats or early termination
if (requireNamespace("matrixStats", quietly = TRUE)) {
  sdnonzero <- matrixStats::colSds(trdat, na.rm = TRUE) > 0
} else {
  # Fallback with colVars approximation
  col_means <- colMeans(trdat, na.rm = TRUE)
  sdnonzero <- colMeans((trdat - rep(col_means, each = nrow(trdat)))^2, na.rm = TRUE) > 0
}
```

---

### 3. `rowwise() + do()` in `extract_prediction_table()`

**Location**: `R/regional.R:31-53, 56-63`

```r
results %>% dplyr::rowwise() %>% dplyr::do({
  # ... complex per-row operations
})
```

**Problem**:
- `do()` is deprecated and slow
- Creates a new tibble for each row then rbinds
- Heavy object allocation overhead

**Impact**: Called for every region in regional analysis.

**Recommended Fix**:
```r
# Use purrr::map_dfr or lapply + bind_rows
purrr::map_dfr(seq_len(nrow(results)), function(i) {
  result <- results$result[[i]]
  id <- results$id[[i]]
  # ... rest of logic
})
```

---

### 4. Repeated Mask Index Computation

**Location**: `R/resample.R:234-237`

```r
if (is.null(mask_indices)) {
  mask_indices <- compute_mask_indices(data$mask)
}
```

**Problem**: Although there's caching at the dataset level (`dset$mask_indices`), if not precomputed, `compute_mask_indices()` does `which(vals != 0)` which is O(n) for n voxels.

**Status**: Partially fixed - the code already precomputes in `mvpa_iterate:399-401`. Ensure all entry points do this.

---

### 5. Cross-validation Sample Generation Inside Inner Loop

**Location**: `R/mvpa_iterate.R:228-229`

```r
samples <- crossval_samples(mspec$crossval,
    tibble::as_tibble(neuroim2::values(roi$train_roi), .name_repair=.name_repair),
    y_train(mspec))
```

**Problem**:
- `crossval_samples()` is called for every ROI
- The CV fold indices are identical across ROIs (same y, same block structure)
- Only the data changes, not the fold assignment

**Impact**: For 50,000 ROIs with 5 folds, this creates 250,000 `modelr::resample` objects.

**Recommended Fix**:
```r
# Precompute fold indices once per analysis
precomputed_folds <- crossval_indices(mspec$crossval, y_train(mspec))

# In internal_crossval, just use the precomputed indices
internal_crossval <- function(mspec, roi, id, center_global_id = NA, fold_indices) {
  data_matrix <- neuroim2::values(roi$train_roi)

  ret <- lapply(seq_along(fold_indices$train), function(i) {
    train_idx <- fold_indices$train[[i]]
    test_idx <- fold_indices$test[[i]]
    train_data <- data_matrix[train_idx, , drop = FALSE]
    test_data <- data_matrix[test_idx, , drop = FALSE]
    # ... rest of training
  })
}
```

---

### 6. Tibble Creation Overhead in Error Paths

**Location**: Multiple locations throughout `mvpa_iterate.R` and `searchlight.R`

```r
tibble::tibble(
  result = list(NULL),
  indices = list(NULL),
  performance = list(NULL),
  id = rnum,
  error = TRUE,
  error_message = "..."
)
```

**Problem**: `tibble()` has more overhead than base R `list()` or `data.frame()`. In error-heavy scenarios (many filtered ROIs), this adds up.

**Recommended Fix**: Use a pre-allocated structure or simple list:
```r
# Create once, reuse
ERROR_TEMPLATE <- list(
  result = list(NULL),
  indices = list(NULL),
  performance = list(NULL),
  id = NA_integer_,
  error = TRUE,
  error_message = ""
)

# In error path
err <- ERROR_TEMPLATE
err$id <- rnum
err$error_message <- msg
```

---

### 7. Sparse Matrix Construction in `combine_randomized`

**Location**: `R/searchlight.R:617-687`

**Current approach**:
```r
for (i in seq_len(nrow(good_results))) {
  # Build triplet lists iteratively
  I_list[[k]] <- rep(ind_i, each = ncols)
  J_list[[k]] <- rep.int(seq_len(ncols), times = len)
  X_list[[k]] <- rep(perf_vec, times = len)
}
```

**Problem**: List appending in a loop can cause repeated memory reallocations.

**Status**: Current implementation is reasonable with pre-allocation via `list()`. Could be marginally improved with `data.table` or Rcpp.

---

## Medium Priority Issues

### 8. `purrr::pmap` vs Base R `Map`

**Location**: `R/mvpa_iterate.R:263-292`

```r
ret <- samples %>% pmap(function(ytrain, ytest, train, test, .id) { ... })
```

**Problem**: `purrr::pmap` has overhead from:
- Type checking and coercion
- Progress tracking infrastructure
- Formula/function parsing

**For hot inner loops**, base R `Map()` or `mapply()` can be 2-3x faster.

---

### 9. Repeated `futile.logger` Calls

**Location**: Throughout all files

```r
futile.logger::flog.debug("message %s", value)
```

**Problem**: Even when log level is INFO, the string formatting still occurs before the level check.

**Recommended Fix**:
```r
# Check level before formatting
if (flog.threshold() <= DEBUG) {
  flog.debug("message %s", value)
}
```

---

## Optimization Priority Matrix

| Issue | Impact | Effort | Priority |
|-------|--------|--------|----------|
| rowwise() in ROI extraction | High | Low | **P0** |
| apply() in filter_roi | High | Low | **P0** |
| CV samples inside loop | High | Medium | **P1** |
| rowwise()+do() in regional | Medium | Low | **P1** |
| Error path tibble creation | Low | Low | **P2** |
| pmap vs Map | Low | Low | **P2** |
| Logger overhead | Low | Low | **P3** |

---

## Benchmarking Suggestions

To measure improvements, use this pattern:

```r
library(bench)

# Benchmark filter_roi improvement
bench::mark(
  original = apply(mat, 2, function(v) any(is.na(v))),
  vectorized = colSums(is.na(mat)) > 0,
  check = TRUE
)

# Profile full searchlight run
profvis::profvis({
  run_searchlight(model_spec, radius = 8, method = "standard")
})
```

---

## Dependencies to Consider

For maximum performance gains:
- `matrixStats` - vectorized column operations
- `data.table` - fast data manipulation
- `Rcpp` - for truly hot inner loops like filter_roi

The `matrixStats` package alone could provide 5-10x speedup for the `filter_roi` operations with minimal code changes.
