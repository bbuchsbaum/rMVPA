# Generate Cross-Validation Folds

Convenience wrapper around
[`crossval_samples`](http://bbuchsbaum.github.io/rMVPA/reference/crossval_samples.md)
for use inside
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md)
methods. Separates fold generation from the iteration engine so that
models can generate folds directly.

## Usage

``` r
generate_folds(cv_spec, data, y)
```

## Arguments

- cv_spec:

  A cross-validation specification object (e.g.,
  `blocked_cross_validation`).

- data:

  A data.frame or tibble of training data.

- y:

  Response variable (factor, numeric vector, or matrix).

## Value

A tibble with columns `ytrain`, `ytest`, `train`, `test`, and `.id`.

## See also

[`crossval_samples`](http://bbuchsbaum.github.io/rMVPA/reference/crossval_samples.md),
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md)
