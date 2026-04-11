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
if (FALSE) { # \dontrun{
# For AUC-0.5 maps (chance-centered)
prepped <- nmf_preprocess_maps(auc_maps, method = "auc")
result <- spatial_nmf_maps(prepped$maps, mask = mask, k = 5)

# For raw AUC maps
prepped <- nmf_preprocess_maps(auc_maps, method = "auc_raw")

# For z-score maps with small positive floor
prepped <- nmf_preprocess_maps(zmaps, method = "shift", min_val = 0.01)
} # }
```
