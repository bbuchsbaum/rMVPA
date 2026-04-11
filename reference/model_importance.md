# Per-Feature Model Importance

Generic function that extracts a per-feature importance vector from a
fitted model object.

## Usage

``` r
model_importance(object, X_train, ...)

# S3 method for class 'sda'
model_importance(object, X_train, summary_fun = NULL, ...)

# S3 method for class 'glmnet'
model_importance(object, X_train, summary_fun = NULL, ...)

# S3 method for class 'spacenet_fit'
model_importance(object, X_train, summary_fun = NULL, ...)

# S3 method for class 'randomForest'
model_importance(object, X_train, ...)

# Default S3 method
model_importance(object, X_train, ...)
```

## Arguments

- object:

  A fitted model object (e.g., from `sda`, `glmnet`, `randomForest`).

- X_train:

  The training data matrix (T x P) used to fit `object`. Required for
  Haufe-based methods; ignored by tree-based methods.

- ...:

  Additional arguments passed to methods.

- summary_fun:

  Optional function to summarize the activation pattern matrix rows.
  Default NULL uses L2 norm.

## Value

A numeric vector of length P with per-feature importance scores, or
`NULL` if importance is not available for this model class.

## Details

The importance measure returned depends on the model class:

- **SDA / glmnet (linear models)**:

  Computes Haufe et al. (2014) *forward-model* activation patterns: A =
  Sigma_x \* W \* inv(W' Sigma_x W). The returned vector is the L2 norm
  (or custom `summary_fun`) of the rows of A. This is the recommended
  importance measure for neuroscience interpretation because it reflects
  where the neural signal originates, not merely which features carry
  discriminative weight.

- **randomForest**:

  Returns MeanDecreaseGini (or MeanDecreaseAccuracy when available).
  This is a **backward** (decoding) measure and is **not** Haufe-safe:
  suppressor variables that reduce node impurity without carrying neural
  signal will receive high importance. Use
  [`region_importance`](http://bbuchsbaum.github.io/rMVPA/reference/region_importance.md)
  for a slower but more robust backward measure that is bounded by
  out-of-sample performance, or restrict neuroscience interpretation to
  linear models with Haufe-based importance.

- **default**:

  Returns `NULL`, signaling that no importance is available for the
  model class.

`model_importance.randomForest` returns Gini-based importance
(MeanDecreaseGini) by default, or permutation-based importance
(MeanDecreaseAccuracy) when the forest was trained with
`importance = TRUE`. Both are **backward** (decoding) measures and may
assign high importance to suppressor variables that do not carry neural
signal. For neuroscience interpretation, prefer
[`haufe_importance`](http://bbuchsbaum.github.io/rMVPA/reference/haufe_importance.md)
with a linear model or use
[`region_importance`](http://bbuchsbaum.github.io/rMVPA/reference/region_importance.md)
for a model-agnostic alternative bounded by out-of-sample performance.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(5,5,5), 40, nlevels=2)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
  vox <- which(ds$dataset$mask > 0)
  X <- neuroim2::series(ds$dataset$train_data, vox)
  fit <- train_model(mspec, X, ds$design$y_train, indices=vox)
  imp <- model_importance(fit, X)
# }
```
