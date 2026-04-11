# Summarize Available Searchlight Engines

Returns a compact table of registered searchlight engines. When a
`model_spec` is supplied, the output also includes whether each engine
is currently eligible for that analysis.

## Usage

``` r
searchlight_engines(
  model_spec = NULL,
  method = c("standard", "randomized", "resampled")
)
```

## Arguments

- model_spec:

  Optional model specification.

- method:

  Searchlight method to audit.

## Value

A data frame.
