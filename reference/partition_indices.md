# Get independent partition indices for a cross-validation fold

Returns the samples belonging to the held-out/independent partition for
a fold. This is distinct from
[`train_indices`](http://bbuchsbaum.github.io/rMVPA/reference/train_indices.md),
which returns the training samples for predictive cross-validation.
Crossnobis estimators use partition indices because their cross-fold
products require independent pattern estimates.

## Usage

``` r
partition_indices(obj, fold_num, ...)
```

## Arguments

- obj:

  A `cross_validation` object.

- fold_num:

  Integer fold number.

- ...:

  Additional arguments passed to methods.

## Value

Integer sample indices for the requested independent partition.
