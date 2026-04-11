# Diagnose Null Distribution Stationarity

Checks whether the null metric values are systematically associated with
covariates, which would indicate that covariate adjustment is needed.

## Usage

``` r
diagnose_null(null_values, covariates, n_perm)
```

## Arguments

- null_values:

  Numeric vector of null metric values.

- covariates:

  A `data.frame` with at least an `nfeatures` column.

- n_perm:

  Integer. Number of permutations used (informational).

## Value

An S3 object of class `"null_diagnostics"`.
