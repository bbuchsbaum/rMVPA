# Pattern-First Spatial Reduced-Rank MVPA Model

Creates a `pattern_model` specification. The model fits a low-rank
forward model \\x = A C^\top y + \epsilon\\ with a structured residual
covariance and derives condition classification, multivariate target
decoding, and encoding predictions from that single fit. It works with
categorical targets (an
[`mvpa_design`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
built from `y_train`) and with vector-valued continuous targets
(`targets = ` a numeric matrix, a
[`feature_sets_design`](https://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md),
or a
[`feature_rsa_design`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md)),
read through
[`model_targets`](https://bbuchsbaum.github.io/rMVPA/reference/model_targets.md).

## Usage

``` r
pattern_model(
  dataset,
  design,
  crossval = NULL,
  rank = "auto",
  max_rank = 8L,
  penalty = NULL,
  noise = list(type = "diag_lowrank", rank = "auto"),
  control = NULL,
  return_predictions = FALSE,
  return_fits = FALSE,
  refit = FALSE,
  ...
)
```

## Arguments

- dataset:

  An
  [`mvpa_dataset`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
  (image, multibasis, surface, or clustered).

- design:

  A design with a
  [`model_targets`](https://bbuchsbaum.github.io/rMVPA/reference/model_targets.md)
  method.

- crossval:

  A cross-validation specification. Defaults to
  [`blocked_cross_validation`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md)
  on the design's block variable, or 5-fold cross-validation when no
  blocks are available.

- rank:

  `"auto"` (selected by nested cross-validation) or a fixed positive
  integer, capped at the eligible rank.

- max_rank:

  Largest rank considered when `rank = "auto"`.

- penalty:

  Spatial penalty specification, a list with elements `sparse`,
  `signed_smooth`, and `support_smooth`. Reserved: this version accepts
  only `NULL`.

- noise:

  Residual covariance specification passed to
  [`pattern_control`](https://bbuchsbaum.github.io/rMVPA/reference/pattern_control.md).

- control:

  A
  [`pattern_control`](https://bbuchsbaum.github.io/rMVPA/reference/pattern_control.md)
  object. When supplied it overrides `max_rank` and `noise`.

- return_predictions:

  Retain out-of-fold prediction ledgers per ROI.

- return_fits:

  Retain the fold fits per ROI (implies ledgers).

- refit:

  Also fit the model on all training rows (a descriptive deployment fit,
  distinct from the cross-validated evidence).

- ...:

  Additional fields stored on the specification.

## Value

A `pattern_model` specification (class
`c("pattern_model", "model_spec")`).

## What is fitted

Targets are centred and whitened on the training rows; \\C\\ has
orthonormal columns so all scale lives in the spatial patterns \\A\\.
The residual covariance \\\Psi = D + UU^\top\\ is estimated once from a
training-only pilot fit and held fixed. Rank is chosen by nested,
block-aware cross-validation on held-out decoding loss (log loss for
categorical targets, normalized squared error for continuous targets).
Spatial penalties are declared through `penalty` but are not yet
available in this version: any non-`NULL` penalty is rejected.

## Outputs

Regional and searchlight runs report scalar metrics per ROI (see
[`output_schema`](https://bbuchsbaum.github.io/rMVPA/reference/output_schema.md)):
`Accuracy`, `AUC`, `logloss`, and `rank_mean` for categorical targets;
`R2`, `RMSE`, `cor`, and `rank_mean` for continuous targets. `rank_mean`
is the mean rank selected across folds and need not be a whole number.
With `return_predictions = TRUE` a regional run also returns the usual
out-of-fold prediction table (categorical targets), and with either
`return_predictions = TRUE` or `return_fits = TRUE` each ROI keeps its
prediction ledger (plus the fold fits when `return_fits = TRUE`) under
`$fits` of the regional result.
[`run_global`](https://bbuchsbaum.github.io/rMVPA/reference/run_global.md)
fits the whole domain once and returns a `pattern_global_result`.

Two ledgers are kept. The fold-resolved one records every prediction
with the fold that produced it. The pooled one holds a single record per
tested observation, in sorted order, with repeats averaged; it is what
the metrics and the prediction table are built from, so a
cross-validation scheme that tests a row several times
([`twofold_blocked_cross_validation`](https://bbuchsbaum.github.io/rMVPA/reference/twofold_blocked_cross_validation.md),
sequential, or bootstrap) does not give that row extra weight. This
matches how every other rMVPA model aggregates repeated predictions.

## See also

[`run_global`](https://bbuchsbaum.github.io/rMVPA/reference/run_global.md),
[`run_regional`](https://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md),
[`predict.pattern_fit`](https://bbuchsbaum.github.io/rMVPA/reference/predict.pattern_fit.md),
[`pattern_control`](https://bbuchsbaum.github.io/rMVPA/reference/pattern_control.md)

## Examples

``` r
ds <- gen_sample_dataset(c(6, 6, 4), 60, nlevels = 3, blocks = 3)
spec <- pattern_model(ds$dataset, ds$design, max_rank = 2)
res <- run_global(spec)
res$performance_table
#> # A tibble: 1 × 4
#>   Accuracy    AUC logloss rank_mean
#>      <dbl>  <dbl>   <dbl>     <dbl>
#> 1      0.3 -0.129    3.32      1.33
```
