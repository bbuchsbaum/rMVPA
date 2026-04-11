# Construct a Global MVPA Result

Construct a Global MVPA Result

## Usage

``` r
global_mvpa_result(
  performance_table,
  result,
  importance_map,
  importance_vector,
  activation_patterns,
  raw_weights,
  fold_fits,
  model_spec
)
```

## Arguments

- performance_table:

  Tibble of cross-validated performance metrics.

- result:

  The merged classification/regression result object.

- importance_map:

  Spatial object (NeuroVol/NeuroSurface) of per-feature importance.

- importance_vector:

  Numeric vector of per-feature importance.

- activation_patterns:

  The P x D activation pattern matrix A.

- raw_weights:

  The averaged P x D weight matrix W.

- fold_fits:

  Optional list of per-fold model_fit objects.

- model_spec:

  The input model specification.

## Value

An S3 object of class `global_mvpa_result`.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Typically created by run_global(), not directly
  result <- global_mvpa_result(
    performance_table = tibble::tibble(Accuracy = 0.8),
    result = NULL,
    importance_map = NULL,
    importance_vector = c(0.1, 0.2),
    activation_patterns = NULL,
    raw_weights = NULL,
    fold_fits = NULL,
    model_spec = list()
  )
} # }
```
