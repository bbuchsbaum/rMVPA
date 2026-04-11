# Construct a Region Importance Result

Construct a Region Importance Result

## Usage

``` r
region_importance_result(
  importance,
  importance_map,
  p_values,
  p_value_map,
  stats_table,
  iteration_log,
  model_spec,
  importance_vector_raw = NULL
)
```

## Arguments

- importance:

  Numeric vector of per-feature importance scores.

- importance_map:

  Spatial object of importance values.

- p_values:

  Numeric vector of Wilcoxon p-values.

- p_value_map:

  Spatial object of -log10(p) values.

- stats_table:

  Tibble with per-feature statistics.

- iteration_log:

  Tibble with per-iteration performance.

- model_spec:

  The input model specification.

- importance_vector_raw:

  Raw importance vector before spatial mapping (default NULL).

## Value

An S3 object of class `region_importance_result`.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Typically created by region_importance(), not directly
  result <- region_importance_result(
    importance = c(0.1, 0.2),
    importance_map = NULL,
    p_values = c(0.05, 0.01),
    p_value_map = NULL,
    stats_table = tibble::tibble(feature_id = 1:2),
    iteration_log = tibble::tibble(iter = 1:10),
    model_spec = list()
  )
} # }
```
