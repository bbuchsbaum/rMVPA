# Build a Covariate-Adjusted Null Distribution

Bins null metric values by a covariate (e.g., number of features per
center) so that per-voxel p-values can be computed relative to a matched
null.

## Usage

``` r
build_adjusted_null(
  null_values,
  covariates,
  n_bins = 5L,
  method = c("adjusted", "global")
)
```

## Arguments

- null_values:

  Numeric vector of null metric values from permutation runs.

- covariates:

  A `data.frame` with an `nfeatures` column (same length as
  `null_values`).

- n_bins:

  Integer \>= 2. Number of quantile bins.

- method:

  Character. `"adjusted"` (default) uses covariate bins; `"global"`
  places all values in one bin.

## Value

An S3 object of class `"adjusted_null"`.
