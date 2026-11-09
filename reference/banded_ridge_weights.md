# Refit a banded-ridge model on all rows and return its weight matrix

Cross-validated accuracies come from \[run_banded_ridge()\]; this
function answers a different question: \*what\* each response is tuned
to. Every response is refit once on all observations, using the
hyperparameters it selected most often across the outer folds, and the
coefficients are returned as one features-by-responses matrix. That
matrix is the input for contrasts over feature columns and for
principal-component summaries of tuning.

## Usage

``` r
banded_ridge_weights(
  model,
  result,
  standardized = TRUE,
  target_batch_size = model$target_batch_size
)
```

## Arguments

- model:

  A \`banded_ridge_model\`.

- result:

  The result of \`run_banded_ridge(model)\`.

- standardized:

  If \`TRUE\` (default), coefficients are per standard deviation of each
  feature, which makes columns with different units comparable. If
  \`FALSE\`, coefficients are in the units of the feature matrix.

- target_batch_size:

  Number of responses read and fit at once.

## Value

A list of class \`banded_ridge_weights\` with

- weights:

  Numeric matrix, features in rows and responses in columns.

- set:

  Feature-band label of each row.

- voxel_index:

  Mask index of each column.

- hyperparameters:

  One row per response: the selected candidate, its \`alpha\` and band
  weights, and in how many outer folds it was chosen.

- standardized:

  The \`standardized\` argument.

- contrast_weights:

  Contrasts-by-responses matrix: the weights projected onto the
  contrasts declared in \`feature_sets(contrasts = )\`, or \`NULL\` when
  none were declared.

If the model has nuisance regressors they are regressed out of features
and responses (with an intercept) before the refit.

## Details

The weights are descriptive. They are estimated from all rows, so they
must not be used to score predictions on those rows.

## See also

\[banded_ridge_model()\], \[run_banded_ridge()\]
