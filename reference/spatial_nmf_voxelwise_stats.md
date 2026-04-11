# Voxelwise Statistics from Spatial NMF Stability

Convenience function that converts stability summaries into voxelwise z-
and p-value maps (based on bootstrap mean/SD).

## Usage

``` r
spatial_nmf_voxelwise_stats(
  x = NULL,
  stability = NULL,
  map_type = NULL,
  mask = NULL,
  dims = NULL,
  mask_indices = NULL,
  full_length = NULL,
  ref_map = NULL
)
```

## Arguments

- x:

  Optional spatial_nmf_maps_result containing stability results.

- stability:

  Optional spatial_nmf_stability result (overrides x\$stability).

- map_type:

  Map type ("volume" or "surface") when providing stability without
  maps.

- mask:

  Mask object for volumetric maps or surface mask (optional for
  surface).

- dims:

  Optional spatial dimensions for volumetric masks.

- mask_indices:

  Optional mask indices used to vectorize maps.

- full_length:

  Optional full map length for surface/vectorized outputs.

- ref_map:

  Optional reference map to copy metadata from.

## Value

A list with \`z\` and \`p\` component maps.

## Examples

``` r
if (FALSE) { # \dontrun{
  # stats <- spatial_nmf_voxelwise_stats(nmf_result, design_matrix)
} # }
```
