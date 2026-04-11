# Test Set Availability

Determine whether the object contains a separate test set.

## Usage

``` r
has_test_set(obj)

# S3 method for class 'mvpa_design'
has_test_set(obj)
```

## Arguments

- obj:

  Object to query.

## Value

Logical indicating if a test set exists.

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
#> external test
has_test_set(ds$design)
#> [1] TRUE
```
