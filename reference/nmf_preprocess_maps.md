# Preprocess Maps for Spatial NMF

Transforms a list of neuroimaging maps to ensure non-negativity for NMF.
NMF requires non-negative input data; this function provides common
transformations for different data types.

## Usage

``` r
nmf_preprocess_maps(
  maps,
  method = c("shift", "auc", "auc_raw", "zscore", "relu", "abs"),
  min_val = 0,
  floor = -0.5,
  mask = NULL
)
```

## Arguments

- maps:

  A list of NeuroVol or NeuroSurface objects.

- method:

  Preprocessing method:

  "shift"

  :   Shifts all values so the minimum becomes \`min_val\` (default 0).
      Use for data with arbitrary negative values.

  "auc"

  :   For AUC values centered at chance (i.e., AUC - 0.5, ranging from
      -0.5 to 0.5). Shifts by 0.5 so chance becomes 0 and perfect
      classification becomes 0.5. Values below \`floor\` (default -0.5)
      are clamped.

  "auc_raw"

  :   For raw AUC values (0 to 1). Subtracts 0.5 then applies "auc"
      method, so chance (0.5) becomes 0.

  "zscore"

  :   For z-scored data. Shifts by \`abs(min) + min_val\`.

  "relu"

  :   Clamps negative values to zero (rectified linear).

  "abs"

  :   Takes absolute value of all data.

- min_val:

  Minimum value after transformation (default 0). A small positive value
  (e.g., 0.01) can help numerical stability.

- floor:

  For "auc" method, values below this are clamped (default -0.5).

- mask:

  Optional mask; if provided, statistics are computed only within mask.

## Value

A list with:

- `maps`: Transformed maps (same class as input).

- `offset`: The offset added (for "shift", "auc", "zscore" methods).

- `method`: The method used.

- `original_range`: Range of original data within mask.

## Details

For group-level NMF analyses, the transformation is computed across all
subjects jointly to preserve relative differences. The returned
\`offset\` can be used to interpret results in the original scale.

## Examples

``` r
set.seed(1)
sp <- neuroim2::NeuroSpace(c(5, 5, 2))
mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(5, 5, 2)), sp)
# Chance-centered AUC maps (AUC - 0.5)
auc_maps <- lapply(1:6, function(i) {
  neuroim2::NeuroVol(array(rnorm(50, sd = 0.1), c(5, 5, 2)), sp)
})

# For AUC-0.5 maps (chance-centered)
prepped <- nmf_preprocess_maps(auc_maps, method = "auc")
result <- spatial_nmf_maps(prepped$maps, mask = mask, k = 2)
#> INFO [2026-10-08 15:55:28] spatial_nmf_maps: fitting NMF (n=6, p=50, k=2, lambda=0)
#> Warning: did not converge--results might be invalid!; try increasing maxit or work
#> INFO [2026-10-08 15:55:29] spatial_nmf_maps: NMF fit complete (converged=TRUE, iterations=47)
#> INFO [2026-10-08 15:55:29] spatial_nmf_maps: parallel=FALSE (explicit=FALSE)

# For raw AUC maps
raw_auc_maps <- lapply(auc_maps, function(m) m + 0.5)
prepped <- nmf_preprocess_maps(raw_auc_maps, method = "auc_raw")

# For z-score maps with small positive floor
zmaps <- lapply(1:6, function(i) neuroim2::NeuroVol(array(rnorm(50), c(5, 5, 2)), sp))
prepped <- nmf_preprocess_maps(zmaps, method = "shift", min_val = 0.01)
prepped$offset
#> [1] 3.018049
```
