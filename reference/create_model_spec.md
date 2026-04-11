# Create a Generic Model Specification

Canonical constructor for custom analysis/model specifications used by
rMVPA's S3 dispatch.

## Usage

``` r
create_model_spec(
  name,
  dataset,
  design,
  return_predictions = FALSE,
  compute_performance = FALSE,
  tune_reps = FALSE,
  ...
)
```

## Arguments

- name:

  A non-empty character scalar indicating the model class name.

- dataset:

  An object inheriting from `mvpa_dataset`.

- design:

  A design object, typically `mvpa_design` or a specialized `*_design`
  class.

- return_predictions:

  Logical scalar indicating whether row-wise predictions should be
  retained.

- compute_performance:

  Logical scalar indicating whether performance metrics should be
  computed during ROI iteration.

- tune_reps:

  Tuning replication control used by model implementations.

- ...:

  Additional named fields stored in the model specification. If
  `has_test_set` is not provided, it is inferred from
  `dataset$test_data` and `design$y_test`.

## Value

A list with class `c(name, "model_spec", "list")`.
