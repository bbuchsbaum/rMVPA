# Describe the frozen target coordinates of a pattern fit

Describe the frozen target coordinates of a pattern fit

## Usage

``` r
pattern_basis(fit)
```

## Arguments

- fit:

  A pattern fit, view, or global result with a refit. Views use the
  underlying fitted basis, not their independent display rotations.

## Value

A `pattern_basis` containing the raw-target-to-score matrix, target
IDs/type, centering, and a hash of the complete frozen transform.

## Details

The intercept in confirmation absorbs target centering. The matrix still
includes training scaling and whitening. Equality of C alone does not
establish compatible coordinates. Target names and physical units must
have the same meaning in every subject.
