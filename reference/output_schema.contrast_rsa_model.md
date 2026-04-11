# Output schema for contrast_rsa_model

Returns a named list describing the output metrics and their lengths.
Each metric is either "scalar" or "vector\[N\]" where N is the number of
contrasts.

## Usage

``` r
# S3 method for class 'contrast_rsa_model'
output_schema(model)
```

## Arguments

- model:

  A `contrast_rsa_model` object.

## Value

A named list where values are "scalar" or "vector\[N\]".
