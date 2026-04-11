# Get Number of Observations

Retrieve the number of observations in an object.

## Usage

``` r
nobs(x)
```

## Arguments

- x:

  The input object.

## Value

The number of observations.

## Examples

``` r
ds <- gen_sample_dataset(c(5,5,5), 20)
nobs(ds$dataset)
#> [1] 20
```
