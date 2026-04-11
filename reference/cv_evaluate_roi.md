# Evaluate One ROI with rMVPA Cross-Validation Helpers

Public wrapper around rMVPA's internal ROI cross-validation helpers.
This allows plugin authors to reuse standard fold execution without
calling non-exported functions.

## Usage

``` r
cv_evaluate_roi(
  model_spec,
  roi_data,
  context = NULL,
  mode = c("auto", "internal", "external"),
  return = c("roi_result", "row"),
  ...
)
```

## Arguments

- model_spec:

  A model specification object.

- roi_data:

  ROI payload as used by
  [`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md).

- context:

  Context list (defaults to
  [`mock_context()`](http://bbuchsbaum.github.io/rMVPA/reference/mock_context.md)
  based on `model_spec`).

- mode:

  One of `"auto"`, `"internal"`, or `"external"`. `"auto"` chooses
  `"external"` when `has_test_set(model_spec)` is true, otherwise
  `"internal"`.

- return:

  One of `"roi_result"` (default) or `"row"`. `"row"` returns the raw
  tibble row produced by the CV helper.

- ...:

  Additional arguments passed to the underlying CV helper.

## Value

A
[`roi_result`](http://bbuchsbaum.github.io/rMVPA/reference/roi_result.md)
by default, or a tibble row when `return = "row"`.
