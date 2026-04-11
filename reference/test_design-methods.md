# Test Design Extraction

Return the design table associated with the test set from an object.

## Usage

``` r
test_design(obj)

# S3 method for class 'mvpa_design'
test_design(obj)
```

## Arguments

- obj:

  The object from which to extract the test design table.

## Value

A data frame containing the test set design variables.

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
#> external test
test_design(ds$design)
#> # A tibble: 10 × 2
#>    Ytest .rownum
#>    <fct>   <int>
#>  1 b           1
#>  2 a           2
#>  3 e           3
#>  4 c           4
#>  5 a           5
#>  6 e           6
#>  7 b           7
#>  8 c           8
#>  9 d           9
#> 10 d          10
```
