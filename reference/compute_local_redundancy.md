# Compute Local Redundancy for Each Searchlight Center

For each searchlight sphere, computes the mean absolute correlation
between the center voxel (first column of the extracted data matrix) and
all other voxels in the sphere. This serves as an autocorrelation proxy
that can be used as a covariate when building the null distribution.

## Usage

``` r
compute_local_redundancy(dataset, searchlight)
```

## Arguments

- dataset:

  An `mvpa_dataset`.

- searchlight:

  A searchlight iterator (list of integer voxel-index vectors).

## Value

Named numeric vector of redundancy values, one per center.
