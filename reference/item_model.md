# ITEM decoding model for ROI/searchlight analysis

Integrates ITEM-style trial-wise decoding into the \`fit_roi\`
architecture by delegating trial-level estimation and covariance-aware
decoding to \`fmrilss\`.

## Usage

``` r
item_model(
  dataset,
  design,
  mode = c("classification", "regression"),
  metric = NULL,
  ridge_u = 0,
  ridge_w = 1e-04,
  lsa_method = c("r", "cpp"),
  solver = c("chol", "svd", "pinv"),
  u_storage = c("matrix", "by_run"),
  class_levels = NULL,
  check_hash = FALSE,
  return_predictions = TRUE,
  compute_performance = TRUE,
  ...
)
```

## Arguments

- dataset:

  An \`mvpa_dataset\`.

- design:

  An \`item_design\` object.

- mode:

  Decoding mode: \`"classification"\` or \`"regression"\`.

- metric:

  Optional ITEM metric. Classification: \`"accuracy"\`,
  \`"balanced_accuracy"\`. Regression: \`"correlation"\`, \`"rmse"\`.

- ridge_u:

  Non-negative ridge used when computing \`U\`.

- ridge_w:

  Non-negative ridge used when fitting ITEM weights per fold.

- lsa_method:

  LS-A backend for \`fmrilss::lsa()\` (\`"r"\` or \`"cpp"\`).

- solver:

  Solver preference for ITEM linear solves (\`"chol"\`, \`"svd"\`,
  \`"pinv"\`).

- u_storage:

  Store trial covariance as full matrix (\`"matrix"\`) or run blocks
  (\`"by_run"\`).

- class_levels:

  Optional class order for classification.

- check_hash:

  Validate trial hash before CV when available.

- return_predictions:

  Keep trial-level prediction tables in ROI results.

- compute_performance:

  Reserved for API compatibility. ITEM always computes scalar ROI
  metrics for mapping.

- ...:

  Additional fields stored on the model spec.

## Value

A model spec of class \`item_model\` compatible with \`run_regional()\`
and \`run_searchlight()\`.

## Details

Conceptually, ITEM differs from \`hrfdecoder_model()\` in where the
decoder is fit: - ITEM (\`item_model\`) first estimates trial-wise
responses (\`Gamma\`) via LS-A and then performs covariance-aware
decoding on those trial estimates. - \`hrfdecoder_model\` fits a
continuous-time decoder directly on TR-level data and aggregates TR
predictions to events afterward.

Use ITEM when you want an explicit trial-estimation stage and direct
control over trial covariance handling (\`U\`), especially for
trial-level diagnostics.
