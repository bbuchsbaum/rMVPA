# Validate a Plugin Model Contract

Executes one
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md) call
and checks that the model returns a valid
[`roi_result`](http://bbuchsbaum.github.io/rMVPA/reference/roi_result.md)
with metrics that match
[`output_schema`](http://bbuchsbaum.github.io/rMVPA/reference/output_schema.md)
(when present).

## Usage

``` r
validate_plugin_model(
  model_spec,
  roi_data = mock_roi_data(),
  context = NULL,
  check_schema = TRUE
)

# S3 method for class 'plugin_validation_result'
print(x, ...)
```

## Arguments

- model_spec:

  A model specification object.

- roi_data:

  ROI payload passed to `fit_roi`. Defaults to
  [`mock_roi_data()`](http://bbuchsbaum.github.io/rMVPA/reference/mock_roi_data.md).

- context:

  Context list passed to `fit_roi`. If `NULL`, a context is constructed
  from `model_spec`.

- check_schema:

  Logical; if `TRUE`, enforce metric name/width agreement with
  `output_schema(model_spec)` when schema is defined.

- x:

  A `plugin_validation_result` object.

- ...:

  Unused.

## Value

`validate_plugin_model()` returns an object of class
`plugin_validation_result`.

`print.plugin_validation_result()` returns `x` invisibly.
