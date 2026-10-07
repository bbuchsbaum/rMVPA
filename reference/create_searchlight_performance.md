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
vol <- neuroim2::NeuroVol(array(runif(27), c(3, 3, 3)),
                          neuroim2::NeuroSpace(c(3, 3, 3)))
perf <- create_searchlight_performance(vol, "Accuracy")
str(perf$summary_stats)
#> List of 4
#>  $ mean: num 0.511
#>  $ sd  : num 0.305
#>  $ min : num 0.0337
#>  $ max : num 0.999
```
