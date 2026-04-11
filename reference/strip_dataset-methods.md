# Strip Dataset from Model Specification

Removes the potentially large dataset component from a model
specification object to avoid copying it during parallel processing.

## Usage

``` r
strip_dataset(obj, ...)

# Default S3 method
strip_dataset(obj, ...)

# S3 method for class 'mvpa_model'
strip_dataset(obj, ...)
```

## Arguments

- obj:

  The model specification object.

- ...:

  Additional arguments.

## Value

The model specification object with the \`dataset\` element removed or
set to NULL.

## Note

For internal parallel dispatch, `as_worker_spec()` is now preferred.
`strip_dataset` remains exported for backward compatibility and external
use.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
  stripped <- strip_dataset(mspec)
  is.null(stripped$dataset)
#> [1] TRUE
# }
```
