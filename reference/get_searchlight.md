# Generate Searchlight Iterator

Generate a searchlight iterator suitable for given data.

## Usage

``` r
get_searchlight(obj, ...)
```

## Arguments

- obj:

  The input dataset object.

- ...:

  Additional arguments to methods.

## Value

A searchlight iterator object.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(5,5,5), 20)
  sl <- get_searchlight(ds$dataset, radius=4)
# }
```
