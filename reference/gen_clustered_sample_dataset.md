# Generate a Synthetic Clustered MVPA Dataset

Creates a synthetic `ClusteredNeuroVec` and associated design for
testing.

## Usage

``` r
gen_clustered_sample_dataset(
  D = c(10, 10, 10),
  nobs = 20,
  K = 5,
  nlevels = 2,
  blocks = 3,
  external_test = FALSE
)
```

## Arguments

- D:

  Spatial dimensions (e.g. `c(10, 10, 10)`).

- nobs:

  Number of observations (time points).

- K:

  Number of clusters.

- nlevels:

  Number of category levels.

- blocks:

  Number of cross-validation blocks.

- external_test:

  Whether to generate an external test set.

## Value

A list with `dataset` (mvpa_clustered_dataset) and `design`
(mvpa_design).

## Examples

``` r
ds <- gen_clustered_sample_dataset(K=10, nobs=20, nlevels=2)
```
