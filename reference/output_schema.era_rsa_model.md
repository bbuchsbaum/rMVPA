# Output schema for era_rsa_model

Declares all scalar metrics emitted by `fit_roi.era_rsa_model` so that
`combine_schema_standard` can pre-allocate output maps with consistent
length across ROIs.

## Usage

``` r
# S3 method for class 'era_rsa_model'
output_schema(model)
```

## Arguments

- model:

  An era_rsa_model object.

## Value

A named character vector mapping metric names to `"scalar"`.
