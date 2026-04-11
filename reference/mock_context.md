# Build a Mock fit_roi Context

Creates the standard `context` list passed to
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md)
methods.

## Usage

``` r
mock_context(design = NULL, cv_spec = NULL, id = 1L, center_global_id = NA_integer_)
```

## Arguments

- design:

  Optional design object. If `NULL`, a small synthetic
  [`mvpa_design`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
  is created.

- cv_spec:

  Optional cross-validation specification.

- id:

  ROI identifier to include in the context.

- center_global_id:

  Optional global center index.

## Value

A list with fields `design`, `cv_spec`, `id`, and `center_global_id`.
