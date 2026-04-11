# Construct a Standardized ROI Result

Creates a uniform return type for per-ROI analysis. Every
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md)
method should return an `roi_result`.

## Usage

``` r
roi_result(
  metrics,
  indices,
  id,
  result = NULL,
  error = FALSE,
  error_message = "~",
  warning = FALSE,
  warning_message = "~"
)
```

## Arguments

- metrics:

  Named numeric vector of performance metrics.

- indices:

  Integer vector of voxel/vertex indices for this ROI.

- id:

  Scalar ROI identifier (center voxel ID or region number).

- result:

  Optional detailed result object (e.g., `classification_result`).

- error:

  Logical; `TRUE` if this ROI failed.

- error_message:

  Character error description, or `"~"` if no error.

- warning:

  Logical; `TRUE` if a warning was raised.

- warning_message:

  Character warning description, or `"~"` if none.

## Value

A list of class `"roi_result"`.

## Examples

``` r
# Successful ROI result
res <- roi_result(
  metrics = c(accuracy = 0.85, AUC = 0.9),
  indices = 1:10,
  id = 42
)
res$metrics
#> accuracy      AUC 
#>     0.85     0.90 

# Error ROI result
err <- roi_result(
  metrics = NULL,
  indices = 1:10,
  id = 42,
  error = TRUE,
  error_message = "Too few features"
)
err$error
#> [1] TRUE
```
