# Test Labels/Response Extraction

Extract the test labels or response variable from an object.

## Usage

``` r
y_test(obj)

# S3 method for class 'mvpa_design'
y_test(obj)

# S3 method for class 'mvpa_model'
y_test(obj)

# S3 method for class 'model_spec'
y_test(obj)
```

## Arguments

- obj:

  The object from which to extract the test response variable.

## Value

The test response variable.

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10, external_test = TRUE)
#> external test
y_test(ds$design)
#>  [1] c a d c e e a d b b
#> Levels: a b c d e
```
