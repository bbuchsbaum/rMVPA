# Validate a Model Specification for Plugin Readiness

Performs structural checks on a model spec and optionally executes a
one-ROI dry run.

## Usage

``` r
validate_model_spec(
  model_spec,
  require_schema = FALSE,
  dry_run = TRUE,
  roi_data = mock_roi_data(),
  context = NULL
)

# S3 method for class 'model_spec_validation_result'
print(x, ...)
```

## Arguments

- model_spec:

  A model specification object.

- require_schema:

  Logical; if `TRUE`, fail when `output_schema(model_spec)` is `NULL`.

- dry_run:

  Logical; if `TRUE`, run
  [`validate_plugin_model`](http://bbuchsbaum.github.io/rMVPA/reference/validate_plugin_model.md)
  on mock ROI/context inputs.

- roi_data:

  ROI payload used for `dry_run`.

- context:

  Context used for `dry_run`.

- x:

  A `model_spec_validation_result` object.

- ...:

  Unused.

## Value

`validate_model_spec()` returns an object of class
`model_spec_validation_result`.

`print.model_spec_validation_result()` returns `x` invisibly.
