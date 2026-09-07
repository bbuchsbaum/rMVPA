# Restrict a Spatial Graph to a Subset of Features

Returns the induced subgraph on the selected feature columns, keeping
the feature-to-column alignment. Use it after dropping invalid or
constant columns from the feature matrix so the graph continues to
describe exactly the columns the estimator sees.

## Usage

``` r
restrict_graph(graph, keep)
```

## Arguments

- graph:

  A `spatial_graph`.

- keep:

  Either a logical vector with one entry per feature column or an
  integer vector of column positions to keep (in the order they should
  appear in the restricted matrix).

## Value

A `spatial_graph` over the kept columns. Its `parent_index` field
records the position of each kept column in the parent graph.

## See also

[`spatial_graph`](https://bbuchsbaum.github.io/rMVPA/reference/spatial_graph.md)

## Examples

``` r
ds <- gen_sample_dataset(c(4, 4, 4), 8)
g <- spatial_graph(ds$dataset)
g2 <- restrict_graph(g, seq_len(10))
g2$n_features
#> [1] 10
```
