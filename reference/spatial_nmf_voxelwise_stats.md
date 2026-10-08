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
set.seed(1)
X <- matrix(runif(20 * 50), 20, 50)
fit <- spatial_nmf(X, k = 2)$fit
stab <- spatial_nmf_stability(X = X, fit = fit, n_boot = 5, seed = 1)
#> INFO [2026-10-08 15:55:49] Spatial NMF stability: 5 bootstrap samples (serial)
#> INFO [2026-10-08 15:55:49] Spatial NMF stability: bootstrap complete.

# Attach bootstrap mean/SD component maps on a small 5 x 5 x 2 grid
sp <- neuroim2::NeuroSpace(c(5, 5, 2))
to_vol <- function(v) neuroim2::NeuroVol(array(v, c(5, 5, 2)), sp)
stab$maps <- list(mean = lapply(1:2, function(i) to_vol(stab$mean[i, ])),
                  sd   = lapply(1:2, function(i) to_vol(stab$sd[i, ])))
stats <- spatial_nmf_voxelwise_stats(stability = stab)
names(stats)
#> [1] "z" "p"
```
