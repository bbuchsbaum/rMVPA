# Create Cross-Validation Folds

Generates a list of row indices for k-fold cross-validation. Can perform
stratified sampling if y is a factor.

## Usage

``` r
create_mvpa_folds(y, k = 5, list = TRUE, seed = NULL)
```

## Arguments

- y:

  A vector, typically the response variable.

- k:

  Integer, the number of folds.

- list:

  Logical, if TRUE, return a list of indices for each fold. If FALSE,
  return a vector of fold assignments for each observation.

- seed:

  Optional integer for reproducible fold creation.

## Value

If \`list=TRUE\`, a list of k integer vectors. If \`list=FALSE\`, an
integer vector of fold assignments.
