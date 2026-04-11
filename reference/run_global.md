# Run Global (Whole-Brain) MVPA Analysis

Train a single classifier on ALL features (parcels or voxels) with
cross-validation, and compute per-feature importance via Haufe et al.
(2014) activation patterns.

## Usage

``` r
run_global(model_spec, ...)

# S3 method for class 'mvpa_model'
run_global(
  model_spec,
  X = NULL,
  summary_fun = NULL,
  return_fits = FALSE,
  aggregation = c("mean", "sum", "maxabs"),
  preflight = c("warn", "error", "off"),
  ...
)
```

## Arguments

- model_spec:

  An `mvpa_model` specification.

- ...:

  Additional arguments (currently unused).

- X:

  Optional pre-computed T x P feature matrix. If NULL, extracted from
  `model_spec$dataset` via `get_feature_matrix`.

- summary_fun:

  Function to summarize activation pattern matrix rows into a scalar
  importance per feature. Default: L2 norm.

- return_fits:

  Logical; if TRUE, store per-fold model fits.

- aggregation:

  How to aggregate multi-basis feature importance (default "mean").

- preflight:

  One of `"warn"` (default), `"error"`, or `"off"`.

## Value

A `global_mvpa_result` object.

A `global_mvpa_result` object.

## Architecture TODO

`run_global` currently uses a dedicated global CV/training pipeline
rather than dispatching through
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md).
This is intentional for now; future cleanup may unify global and ROI
fitting interfaces.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(5,5,5), 40, nlevels=2, blocks=3)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
    "classification", crossval=cval)
  result <- run_global(mspec)
# }
```
