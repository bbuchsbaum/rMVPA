# Create an MVPA Dataset for Clustered/Parcellated Data

Creates a dataset object for MVPA analysis on parcellated brain data
stored as `ClusteredNeuroVec` objects from neuroim2. A "clustered
searchlight" iterates over all K parcels, pooling nearby parcels as
features for each analysis.

## Usage

``` r
mvpa_clustered_dataset(train_data, test_data = NULL, mask = NULL)
```

## Arguments

- train_data:

  A `ClusteredNeuroVec` instance (training data).

- test_data:

  An optional `ClusteredNeuroVec` instance (test data). Default: NULL.

- mask:

  An optional mask. If NULL, derived from `train_data`'s cluster volume.

## Value

An `mvpa_clustered_dataset` object (S3 class) containing:

- train_data:

  The training data as a `ClusteredNeuroVec`.

- test_data:

  The test data (if provided).

- mask:

  The mask.

- cvol:

  The `ClusteredNeuroVol` defining cluster assignments.

- region_mask:

  A `NeuroVol` where each voxel's value is its cluster ID (for
  broadcasting results via `map_values`).

- has_test_set:

  Logical flag.

## Examples

``` r
if (FALSE) { # \dontrun{
  ds <- gen_clustered_sample_dataset(K=5, nobs=20)
  print(ds$dataset)
} # }
```
