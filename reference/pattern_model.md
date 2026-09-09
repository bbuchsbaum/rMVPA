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
  graph = NULL,
  noise = list(type = "diag_lowrank", rank = "auto"),
  control = NULL,
  return_predictions = FALSE,
  return_fits = FALSE,
  refit = FALSE,
  weights = NULL,
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

  Spatial penalty on the forward patterns: a list with any of `sparse`
  (row-wise group lasso, as a fraction of the penalty that empties the
  model), `signed_smooth` (graph-Laplacian smoothing of the signed
  loadings, relative to the data-fit curvature), and `support_smooth`
  (reserved, not implemented). Each may be a number, a vector of
  candidates, or `"auto"`. `NULL` (default) fits without a spatial
  penalty. See the Spatial penalties section.

- graph:

  A
  [`spatial_graph`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_graph.md)
  aligned to the dataset's features, used by `signed_smooth`. Built from
  the dataset when needed and not supplied.

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

- weights:

  Optional non-negative observation weights, one per training
  observation. When `NULL` (default) they are read from the design
  ([`feature_sets_design`](https://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md)
  carries `row_weights`); a vector supplied here overrides the design's
  weights. Weights enter the estimator itself – centring, target
  whitening, the residual-covariance estimate, and the penalized
  objective are all weighted – and the held-out loss that drives
  rank/penalty selection. Integer weights are exactly equivalent to
  replicating rows, and a zero weight is exactly equivalent to omitting
  the row from training. Reported performance metrics remain unweighted:
  every tested observation counts once.

- ...:

  Additional fields stored on the specification.

## Value

A `pattern_model` specification (class
`c("pattern_model", "model_spec")`).

## What is fitted

Targets are centred and whitened on the training rows; \\C\\ has
orthonormal columns so all scale lives in the spatial patterns \\A\\.
The residual covariance \\\Psi = D + UU^\top\\ is estimated once from a
training-only pilot fit and held fixed. Rank and penalty strength are
chosen together by nested, block-aware cross-validation on held-out
decoding loss (log loss for categorical targets, normalized squared
error for continuous targets), with ties going to the smaller rank.

## Spatial penalties

Penalties act on the forward patterns \\A\\, not on decoding weights, so
what is regularized is where task signal is expressed rather than which
measurements happen to help prediction.

- `sparse`:

  Row-wise group lasso: a feature is either in the patterns or out of
  them, for all components at once. Given as a fraction of the smallest
  penalty that empties the model, so the same number means the same
  thing across folds and feature domains. Needs no anatomy.

- `signed_smooth`:

  A graph-Laplacian quadratic that pulls neighbouring loadings together,
  weighted relative to the curvature of the data-fit term, so 1 makes
  smoothing as influential as the data. Requires a
  [`spatial_graph`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_graph.md),
  built from the dataset unless one is supplied. This assumes
  neighbouring voxels carry *similar signed* loadings, which is wrong
  for a fine-grained code that flips sign within a region, so it
  defaults to off and is worth tuning.

- `support_smooth`:

  Reserved for a spatially smooth support envelope, which would let a
  coherent anatomical territory contain sign-flipping loadings. Not
  implemented; it is rejected rather than quietly redirected to
  `signed_smooth`, because they encode different assumptions.

Either penalty may be `"auto"`, which cross-validates over a short path;
a number is used as given and does not enlarge the tuning grid.

## Outputs

Regional and searchlight runs report scalar metrics per ROI (see
[`output_schema`](https://bbuchsbaum.github.io/rMVPA/reference/output_schema.md)):
`Accuracy`, `AUC`, `logloss`, and `rank_mean` for categorical targets;
`R2`, `RMSE`, `cor`, and `rank_mean` for continuous targets. `rank_mean`
is the mean rank selected across folds and need not be a whole number.
When a penalty is in force, `n_selected` reports the mean number of
features with a non-zero pattern. With `return_predictions = TRUE` a
regional run also returns the usual out-of-fold prediction table
(categorical targets), and with either `return_predictions = TRUE` or
`return_fits = TRUE` each ROI keeps its prediction ledger (plus the fold
fits when `return_fits = TRUE`) under `$fits` of the regional result.
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

# sparse forward patterns, with the penalty chosen by nested CV
sparse_spec <- pattern_model(ds$dataset, ds$design, rank = 1,
                             penalty = list(sparse = "auto"))
run_global(sparse_spec)$performance_table
#> # A tibble: 1 × 5
#>   Accuracy    AUC logloss rank_mean n_selected
#>      <dbl>  <dbl>   <dbl>     <dbl>      <dbl>
#> 1     0.25 -0.128    1.33         1       13.3
```
