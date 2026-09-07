# Edge List of a Spatial Graph

Edge List of a Spatial Graph

## Usage

``` r
graph_edges(graph)
```

## Arguments

- graph:

  A `spatial_graph`.

## Value

A data frame with integer columns `from` and `to` (column positions,
`from < to`) and a numeric `weight`.

## See also

[`spatial_graph`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_graph.md)

## Examples

``` r
ds <- gen_sample_dataset(c(4, 4, 4), 8)
head(graph_edges(spatial_graph(ds$dataset)))
#>   from to weight
#> 1    1  2      1
#> 2    2  3      1
#> 3    3  4      1
#> 4    1  5      1
#> 5    2  6      1
#> 6    5  6      1
```
