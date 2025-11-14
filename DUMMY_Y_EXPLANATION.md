# Understanding the Dummy `y` Variable in hrfdecoder_design

## Question
Why does `hrfdecoder_design` use `y = seq_len(T)` (TR indices 1:T) as a dummy variable? Is this correct?

## Answer: Yes, it's correct!

The `y = 1:T` pattern is **architecturally sound** for continuous decoding in rMVPA. Here's why:

## How rMVPA's Cross-Validation Works

### Fold Construction (blocked_cross_validation)
```r
# From R/crossval.R, crossv_block():
crossv_block <- function(data, y, block_var, id = ".id") {
  idx <- seq_len(nrow(data))
  fold_idx <- split(idx, block_var)  # ← Splits by block_var, NOT y!

  fold <- function(test) {
    tidx <- setdiff(idx, test)
    list(
      ytrain = subset_y(y, tidx),  # ← Just subsets y by indices
      ytest = subset_y(y, test),
      train = modelr::resample(data, tidx),
      test = modelr::resample(data, test)
    )
  }
}
```

**Critical insight:** Fold assignment is determined **entirely by `block_var`**, not by `y` values!

### What `y` is Actually Used For

1. **Length validation**: `check_len(y, block_var)` ensures `length(y) == length(block_var)`
2. **Subsetting**: Creates `ytrain` and `ytest` by indexing: `y[train_indices]`
3. **Passed to train_model**: But `train_model.hrfdecoder_model()` **explicitly ignores it**

```r
# From R/hrfdecoder_model.R:
train_model.hrfdecoder_model <- function(obj, train_dat, y, ...) {
  # y is IGNORED - not used anywhere in this function!

  # Actual targets come from design metadata:
  evm <- obj$design$event_model
  events <- obj$design$events

  fit <- hrfdecoder::hrfdecoder_fit(
    X = as.matrix(train_dat),
    event_model = evm,  # ← Real temporal structure
    ...
  )
}
```

## Why Use `1:T` Instead of Alternatives?

| Option | Code | Pros | Cons |
|--------|------|------|------|
| **Sequential (current)** | `seq_len(T)` | Explicit, shows T observations | None! |
| NULL | `NULL` | Conceptually clean | Might confuse readers |
| Constant | `rep(1, T)` | Simple | Obscures granularity |

**Choice: `seq_len(T)` is best because:**
- ✅ **Explicit**: Makes TR-level granularity clear
- ✅ **Safe**: Guaranteed to pass length validation
- ✅ **Informative**: Shows there are T observations (TRs)
- ✅ **Conventional**: Matches how indices are typically represented

## The Data Flow

```
User specifies:
├─ event_model (fmridesign): HRF-convolved event design
├─ events (data.frame): Event-level labels (onset, condition)
└─ block_var (vector): Run IDs per TR

hrfdecoder_design creates:
├─ y_train = 1:T  ← DUMMY (for CV machinery only)
├─ block_var = run IDs  ← REAL (determines folds)
└─ Stores event_model + events  ← REAL (decoding targets)

Cross-validation:
├─ Splits folds by block_var (e.g., hold out run 3)
├─ Subsets y → ytrain, ytest (passed but ignored)
└─ Trains on TR-level data using event_model/events

Training:
├─ train_model receives ytrain (ignored)
├─ Uses event_model + events for actual targets
└─ Fits on continuous TRs, aggregates to events
```

## Comparison with Standard MVPA

### Standard MVPA (Trial-Level)
```r
# y_train contains REAL labels (per trial)
y_train <- c("face", "house", "face", "house", ...)  # Length = # trials

mvpa_design(
  y_train = ~ condition,  # Real trial labels
  block_var = ~ run       # Run IDs per trial
)
```
- Folds split by trials (some trials in train, others in test)
- `y_train` is used for training the classifier
- One observation = one trial

### hrfdecoder (Continuous TR-Level)
```r
# y_train is DUMMY (TR indices)
y_train <- 1:200  # Length = # TRs, not # events

hrfdecoder_design(
  event_model = evmod,    # Real event structure
  events = events_df,     # Real event labels
  block_var = run_ids     # Run IDs per TR
)
```
- Folds split by runs (entire runs held out)
- Dummy `y` ignored; real targets in `event_model`/`events`
- One observation = one TR (not one trial!)

## Key Architectural Points

1. **Blocked CV doesn't use `y` for fold assignment**
   - Only `block_var` determines which observations go in which fold
   - `y` is just metadata that gets subset and passed through

2. **hrfdecoder trains on continuous TRs**
   - Operates on all TRs from training runs
   - Real targets (events) come from design metadata
   - Aggregates TR-level predictions to event-level post-hoc

3. **The pattern is well-established**
   - Other specialized models (RSA, MANOVA) also use `y` in non-standard ways
   - The S3 generic system allows models to interpret `y` however they need

## Documentation Improvements Made

### 1. Enhanced inline comments in `hrfdecoder_design.R`:
```r
# IMPORTANT: y is a DUMMY variable (TR sequence 1:T) that serves THREE purposes:
# 1. Length validation: crossval_samples() checks length(y) == nrow(data)
# 2. Subsetting: Creates ytrain/ytest by indexing into y
# 3. Passed to train_model: But train_model.hrfdecoder_model() IGNORES it
#
# The ACTUAL decoding targets come from $event_model and $events.
# Fold assignment is determined SOLELY by block_var (run IDs), NOT by y values.
```

### 2. Updated roxygen2 documentation:
- Clarified that fold assignment uses `block_var` only
- Explained that `train_model` ignores the `y` parameter
- Added explicit "IGNORED" annotation in parameter docs

### 3. Enhanced vignette explanation:
- Clear bullet points explaining dummy `y` purpose
- Emphasized fold assignment mechanism
- Contrasted with standard MVPA approach

## Conclusion

The `y = seq_len(T)` pattern is:
- ✅ **Architecturally correct**
- ✅ **Well-documented** (after improvements)
- ✅ **Following rMVPA conventions**
- ✅ **Making TR-level granularity explicit**

The design is **kosher** - it leverages rMVPA's flexible S3 system exactly as intended!
