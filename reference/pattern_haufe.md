# Compare empirical held-out Haufe patterns with the forward model

Compare empirical held-out Haufe patterns with the forward model

## Usage

``` r
pattern_haufe(fit, X_holdout, observation_ids = NULL)
```

## Arguments

- fit:

  A pattern fit or view.

- X_holdout:

  Independent rows, with all original input columns.

- observation_ids:

  Held-out row identifiers in the same namespace as
  `fit$training_observation_ids`. Required when the fit carries those
  IDs.

## Value

Empirical and model patterns on retained features in original units,
relative Frobenius discrepancy, effective score rank, and row
identifiers.

## Details

Independence is a caller obligation; an ID overlap is rejected. For a
rank-deficient decoder, the comparison uses the identifiable score
subspace `A G+ G`. This is a covariance diagnostic, not an inference
test.
