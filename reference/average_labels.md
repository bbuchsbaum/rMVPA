# Average NeuroVec Data by Labels

Computes the average brain activation pattern for each unique
label/condition, with optional normalization of individual volumes
before averaging.

## Usage

``` r
average_labels(
  neurovec,
  labels,
  mask = NULL,
  normalize = c("none", "scale", "center", "z", "unit", "percent"),
  normalize_by = c("volume", "voxel"),
  na.rm = TRUE,
  return_matrix = FALSE
)
```

## Arguments

- neurovec:

  A NeuroVec object containing the neuroimaging data

- labels:

  A vector of labels/conditions corresponding to each volume in neurovec

- mask:

  An optional mask (NeuroVol or logical array) to restrict analysis to
  specific voxels

- normalize:

  Character string specifying normalization method applied to each
  volume before averaging. Note that after averaging normalized volumes,
  the resulting averages may not maintain the same normalization
  properties (e.g., averaged z-scored data will have SD \< 1):

  "none"

  :   No normalization (default)

  "scale"

  :   Scale each volume to unit variance (divide by SD)

  "center"

  :   Center each volume (subtract mean)

  "z"

  :   Z-score normalization (center and scale)

  "unit"

  :   Scale to unit norm (L2 normalization)

  "percent"

  :   Convert to percent signal change from volume mean

- normalize_by:

  Character string specifying normalization scope:

  "volume"

  :   Normalize within each volume (default)

  "voxel"

  :   Normalize each voxel across time

- na.rm:

  Logical, whether to remove NA values during averaging (default TRUE)

- return_matrix:

  Logical, if TRUE returns the averaged data matrix instead of NeuroVec
  (default FALSE)

## Value

A NeuroVec object with averaged data (one volume per unique label), or a
matrix if return_matrix=TRUE. The returned object has attributes:

- condition_labels:

  The unique labels in order

- n_averaged:

  Number of volumes averaged for each label

- normalization:

  The normalization method used

## Examples

``` r
# Small synthetic NeuroVec: 4x4x4 volume, 12 time points, 3 conditions
sp <- neuroim2::NeuroSpace(c(4, 4, 4, 12))
scandat <- neuroim2::NeuroVec(array(rnorm(4 * 4 * 4 * 12), c(4, 4, 4, 12)), sp)
mask <- neuroim2::LogicalNeuroVol(array(TRUE, c(4, 4, 4)),
                                  neuroim2::NeuroSpace(c(4, 4, 4)))
condition_labels <- rep(c("A", "B", "C"), times = 4)

# Basic averaging
averaged <- average_labels(scandat, condition_labels, mask)

# With z-score normalization of each volume
averaged_norm <- average_labels(scandat, condition_labels, mask,
                                normalize = "z", normalize_by = "volume")

# Scale to unit norm for RSA
averaged_unit <- average_labels(scandat, condition_labels, mask,
                                normalize = "unit")

# Get just the data matrix (conditions x voxels)
data_mat <- average_labels(scandat, condition_labels, mask,
                           return_matrix = TRUE)
dim(data_mat)
#> [1]  3 64
```
