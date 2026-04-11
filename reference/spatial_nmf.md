# Spatial Non-negative Matrix Factorization

A single high-level interface for spatial NMF. Use `spatial_nmf()` with
a non-negative subject-by-voxel matrix, or with a list of
NeuroVol/NeuroSurface maps. Optional preprocessing can make signed maps
or matrices non-negative before factorization.

## Usage

``` r
spatial_nmf(
  x,
  group_B = NULL,
  groups = NULL,
  k,
  mask = NULL,
  dims = NULL,
  transform = c("none", "shift", "auc", "auc_raw", "zscore", "relu", "abs"),
  min_val = 0,
  floor = -0.5,
  lambda = 0,
  fast = FALSE,
  graph = NULL,
  neighbors = 6,
  na_action = c("zero", "error"),
  return_maps = .is_map_input(x),
  return_data = FALSE,
  component_test = NULL,
  global_test = NULL,
  stability = NULL,
  voxelwise_stats = c("none", "stability_zp"),
  parallel = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- x:

  A numeric subject-by-voxel matrix/data frame, a single map, or a list
  of NeuroVol/NeuroSurface maps for group A.

- group_B:

  Optional list of NeuroVol/NeuroSurface maps for group B. Use this only
  when `x` is map input.

- groups:

  Optional group labels for matrix input. Required for
  `component_test = TRUE` or `global_test = TRUE` with matrix input.

- k:

  Number of components.

- mask:

  Mask object for volumetric/surface maps, or optional spatial metadata
  for matrix input when `return_maps = TRUE` or `lambda > 0`.

- dims:

  Optional spatial dimensions for volumetric masks.

- transform:

  Input transformation before NMF: "none" requires non-negative input;
  other values use the same semantics as
  [`nmf_preprocess_maps`](http://bbuchsbaum.github.io/rMVPA/reference/nmf_preprocess_maps.md).

- min_val:

  Minimum value after transformation.

- floor:

  Lower floor used by the "auc" and "auc_raw" transformations.

- lambda:

  Spatial regularization strength (0 = none).

- fast:

  Logical; use faster, lower-iteration defaults in the NMF fit.

- graph:

  Optional graph Laplacian list from `build_graph_laplacian`.

- neighbors:

  Neighborhood size for volumetric adjacency (6/18/26).

- na_action:

  How to handle NA values: "zero" or "error".

- return_maps:

  Logical; return component maps when spatial metadata is available.

- return_data:

  Logical; include the subject-by-voxel data matrix.

- component_test:

  NULL to skip, TRUE for defaults, or a list of arguments passed to
  [`spatial_nmf_component_test`](http://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_component_test.md).

- global_test:

  NULL to skip, TRUE for defaults, or a list of arguments passed to
  [`spatial_nmf_global_test`](http://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_global_test.md).

- stability:

  NULL to skip, TRUE for defaults, or a list of arguments passed to
  [`spatial_nmf_stability`](http://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_stability.md).

- voxelwise_stats:

  Optional voxelwise statistics. Use "stability_zp" to return bootstrap
  z- and p-value maps derived from stability.

- parallel:

  Logical; enable parallel processing for inference functions.

- progress:

  Logical; report progress via progressr package.

- ...:

  Additional arguments passed to the NMF optimizer.

## Value

A `spatial_nmf_result` object with `fit`, optional `components`,
optional `data`, and requested inference results.

## Examples

``` r
set.seed(1)
W <- matrix(runif(12 * 2), nrow = 12)
H <- matrix(runif(2 * 20), nrow = 2)
X <- W %*% H
result <- spatial_nmf(X, k = 2, max_iter = 20)

if (FALSE) { # \dontrun{
result <- spatial_nmf(group_A_vols, group_B = group_B_vols, mask = mask, k = 5)
} # }
```
