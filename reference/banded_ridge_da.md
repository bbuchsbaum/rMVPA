# Convenience wrapper: build a grouped-ridge domain-adaptation model from matrices

\`banded_ridge_da()\` is a convenience wrapper that builds:

- a \`feature_sets\` object for train predictors,

- a test \`feature_sets\` object (from \`X_test\` or from \`gamma\` via
  \`expected_features()\`),

- a \`feature_sets_design\`, and

- the final \`banded_ridge_da_model\` spec.

Preferred name for \`banded_ridge_da()\`. See that function for full
details.

## Usage

``` r
banded_ridge_da(
  dataset,
  X_train,
  spec = NULL,
  X_test = NULL,
  gamma = NULL,
  target_builder = NULL,
  target_builder_data = NULL,
  n_test = NULL,
  drop_null = TRUE,
  renormalize = FALSE,
  block_var_test = NULL,
  ...
)

grouped_ridge_da(
  dataset,
  X_train,
  spec = NULL,
  X_test = NULL,
  gamma = NULL,
  target_builder = NULL,
  target_builder_data = NULL,
  n_test = NULL,
  drop_null = TRUE,
  renormalize = FALSE,
  block_var_test = NULL,
  ...
)
```

## Arguments

- dataset:

  mvpa_dataset with train_data/test_data.

- X_train:

  Train predictor matrix (T_train x D) or a \`feature_sets\` object.

- spec:

  Feature-set spec for matrix inputs, created by \`blocks()\` or
  \`by_set()\`. Ignored if \`X_train\` is already a \`feature_sets\`.

- X_test:

  Optional test predictor matrix (T_test x D) or a \`feature_sets\`
  object.

- gamma:

  Optional alignment matrix used when \`X_test\` is NULL. See
  \`expected_features()\`.

- target_builder:

  Optional fold-aware callback passed through to
  \`feature_sets_design()\`. It can rebuild target predictors separately
  for each outer target fold and may return a \`feature_sets\` object,
  numeric matrix, or a list containing \`gamma\`, \`X\`, or \`X_test\`.

- target_builder_data:

  Optional object passed through to \`target_builder\` as
  \`builder_data\`.

- n_test:

  Optional target row count used when \`target_builder\` is provided
  without a fixed \`X_test\`.

- drop_null, renormalize:

  Passed to \`expected_features()\` when using \`gamma\`.

- block_var_test:

  Optional test run/block vector (length T_test).

- ...:

  Passed through to \`banded_ridge_da_model()\`.

## Value

A model spec of class \`banded_ridge_da_model\`.

A model spec of class \`banded_ridge_da_model\`.

## Details

Use this when you already have `X_train` (TR x features) as a single
matrix and you want to declare sets via \`blocks()\` or \`by_set()\`.

## See also

[`banded_ridge_da_model`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md),
[`grouped_ridge_da_model`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md),
[`blocks`](http://bbuchsbaum.github.io/rMVPA/reference/blocks.md),
[`by_set`](http://bbuchsbaum.github.io/rMVPA/reference/by_set.md)

## Examples

``` r
if (FALSE) { # \dontrun{
ms <- grouped_ridge_da(
  dataset = dset,
  X_train = X_enc,
  spec = blocks(low = 100, mid = 100, high = 100, sem = 100),
  gamma = gamma,
  block_var_test = recall_runs,
  mode = "stacked",
  lambdas = c(low = 10, mid = 10, high = 10, sem = 10),
  alpha_recall = 0.2
)
} # }
```
