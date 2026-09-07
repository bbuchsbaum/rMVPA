# Spatial Adjacency Graph Aligned to Dataset Features

Build a sparse adjacency graph whose vertex \\j\\ refers to column \\j\\
of the feature matrix returned by
[`get_feature_matrix`](https://bbuchsbaum.github.io/rMVPA/reference/get_feature_matrix.md)
for the same dataset. This alignment is the contract every spatially
regularized estimator in rMVPA relies on.

## Usage

``` r
spatial_graph(x, ...)

# S3 method for class 'mvpa_image_dataset'
spatial_graph(x, neighbors = 6, ...)

# S3 method for class 'mvpa_multibasis_image_dataset'
spatial_graph(x, neighbors = 6, connect_basis = FALSE, ...)

# S3 method for class 'mvpa_surface_dataset'
spatial_graph(x, ...)

# S3 method for class 'mvpa_clustered_dataset'
spatial_graph(x, neighbors = 6, ...)

# Default S3 method
spatial_graph(
  x,
  feature_ids = NULL,
  weighted = FALSE,
  domain_type = "custom",
  ...
)
```

## Arguments

- x:

  A dataset
  ([`mvpa_dataset`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
  and friends) or a raw adjacency specification.

- ...:

  Additional arguments passed to methods; see
  `spatial_graph.mvpa_image_dataset`.

- neighbors:

  Voxel neighbourhood for volumetric grids: 6 (faces), 18 (faces and
  edges), or 26 (faces, edges, and corners).

- connect_basis:

  Logical; for multibasis datasets, also connect each voxel to the same
  voxel in every other basis channel (default `FALSE`: channels form
  disconnected components, so smoothing never crosses channels).

- feature_ids:

  Integer identifiers mapping graph vertices to dataset locations when a
  raw adjacency is supplied. Defaults to `seq_len(n)`.

- weighted:

  Logical; keep edge weights of a raw adjacency instead of binarizing
  them.

- domain_type:

  Label for a raw adjacency's domain (default `"custom"`).

## Value

A `spatial_graph` object.

## Details

Methods exist for volumetric datasets (grid adjacency over the mask; 6,
18, or 26 neighbours), multibasis volumetric datasets (one grid graph
per basis channel, disconnected across channels unless
`connect_basis = TRUE`), surface datasets (mesh adjacency from
neurosurf, restricted to masked nodes), clustered datasets (two parcels
are adjacent when any of their voxels are), and raw adjacency matrices
supplied directly (`matrix`, `Matrix`, or a list with an `A` element).

The result is a list of class `spatial_graph` with the sparse symmetric
adjacency `A`, its `degree` vector, the combinatorial Laplacian
`L = diag(degree) - A`, the `feature_ids` that map graph vertices back
to dataset locations, the number of features `n_features`, a
`domain_type`, a `geometry_id` string identifying the geometry the graph
was built from, and, for multibasis data, a `basis` vector giving each
column's basis channel. The object also carries the `A`, `weighted`, and
`degree` fields that the graph-regularized NMF functions accept as their
`graph` argument.

For clustered datasets the graph vertices are the parcels, in the column
order of the cluster time series, and `feature_ids` holds the actual
cluster labels (which need not be `1..K`). Column positions, not labels,
index the feature matrix.

## See also

[`restrict_graph`](https://bbuchsbaum.github.io/rMVPA/reference/restrict_graph.md),
[`graph_edges`](https://bbuchsbaum.github.io/rMVPA/reference/graph_edges.md)

## Examples

``` r
ds <- gen_sample_dataset(c(5, 5, 5), 10)
g <- spatial_graph(ds$dataset)
g$n_features == ncol(get_feature_matrix(ds$dataset))
#> [1] TRUE
```
