# Create Searchlight Performance Object

Creates a searchlight_performance object with the expected structure for
tests

## Usage

``` r
create_searchlight_performance(data, metric_name, indices = NULL)
```

## Arguments

- data:

  NeuroVol or NeuroSurface object

- metric_name:

  Character string naming the metric

- indices:

  Numeric vector of center indices (optional)

## Value

A searchlight_performance object

## Examples

``` r
if (FALSE) { # \dontrun{
  # Internal function for creating searchlight performance objects
  perf <- create_searchlight_performance(neurovol, "Accuracy")
} # }
```
