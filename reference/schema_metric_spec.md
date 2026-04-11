# Metric Schema Constructors

Helpers for defining typed metric entries in
[`output_schema`](http://bbuchsbaum.github.io/rMVPA/reference/output_schema.md)
methods for plugin models.

## Usage

``` r
schema_scalar()

schema_vector(n)
```

## Arguments

- n:

  Positive integer length for vector-valued metrics.

## Value

An object of class `rmvpa_metric_spec`.

## Details

`schema_scalar()` declares a single metric column. `schema_vector(n)`
declares `n` metric columns expanded as `metric.1`, `metric.2`, ...,
`metric.n`.

These constructors are optional. Legacy string specs (`"scalar"`,
`"vector[N]"`) remain supported.

## Examples

``` r
schema_scalar()
#> $type
#> [1] "scalar"
#> 
#> $size
#> [1] 1
#> 
#> attr(,"class")
#> [1] "rmvpa_metric_spec"
schema_vector(3)
#> $type
#> [1] "vector"
#> 
#> $size
#> [1] 3
#> 
#> attr(,"class")
#> [1] "rmvpa_metric_spec"
```
