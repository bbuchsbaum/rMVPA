# Extract forward patterns, calibrated weights, or invariant maps

Extract forward patterns, calibrated weights, or invariant maps

## Usage

``` r
model_patterns(
  fit,
  type = c("forward", "weights", "conditional_info", "signal_sd"),
  dataset = NULL,
  ...
)

# S3 method for class 'pattern_fit'
model_importance(
  object,
  X_train = NULL,
  type = c("signal_sd", "conditional_info"),
  ...
)

# S3 method for class 'pattern_view'
model_importance(
  object,
  X_train = NULL,
  type = c("signal_sd", "conditional_info"),
  ...
)

# S3 method for class 'pattern_global_result'
model_importance(
  object,
  X_train = NULL,
  type = c("signal_sd", "conditional_info"),
  ...
)
```

## Arguments

- fit:

  A `pattern_fit`, `pattern_view`, or global result with a refit.

- type:

  Quantity to extract.

- dataset:

  Optional dataset for `build_output_map`. Global results supply their
  dataset automatically. Multibasis image maps are refused because
  aggregation across channels changes the estimand; extract vectors.

- ...:

  Reserved.

- object:

  A fit, view, or global result.

- X_train:

  Unused; maps are implied by the fitted model.

## Value

A vector or matrix aligned to all input columns (screened columns are
`NA`), or an image / list of component images when a dataset is
supplied.

## Details

Forward patterns and signal standard deviations are in original feature
units. Weights map original, centred measurements to calibrated scores.
Forward patterns and weights depend on the component basis. Signal SD
and conditional information are coordinate invariant. Conditional
information is in nats under the working Gaussian model for task scores,
not empirical information about class labels. A feature with zero
loading can have positive information through noise cancellation.
