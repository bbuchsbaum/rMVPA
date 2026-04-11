# Get Searchlight Scope

Returns the analysis scope for a dataset, controlling how `mvpa_iterate`
handles center IDs and whether `pobserved` maps are constructed.

## Usage

``` r
searchlight_scope(dataset, ...)
```

## Arguments

- dataset:

  The dataset object.

- ...:

  Additional arguments.

## Value

Character string: `"searchlight"` or `"regional"`.
