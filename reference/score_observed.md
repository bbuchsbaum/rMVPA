# Compute Permutation P-Values

For each observed metric value, assigns it to the matching null bin and
computes a conservative one-sided p-value.

## Usage

``` r
score_observed(observed_values, adjusted_null, covariates_full)
```

## Arguments

- observed_values:

  Numeric vector of observed metric values (one per searchlight center
  across the full brain).

- adjusted_null:

  An `"adjusted_null"` object from
  [`build_adjusted_null`](http://bbuchsbaum.github.io/rMVPA/reference/build_adjusted_null.md).

- covariates_full:

  A `data.frame` with an `nfeatures` column for all voxels (same length
  as `observed_values`).

## Value

Numeric vector of p-values in \\\[0, 1\]\\.
