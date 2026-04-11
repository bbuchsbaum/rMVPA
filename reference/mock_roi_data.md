# Build Mock ROI Data for Plugin Tests

Constructs a minimal ROI payload compatible with
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md).
This is useful for unit-testing plugin methods without running
[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
or
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md).

## Usage

``` r
mock_roi_data(
  train_data = NULL,
  test_data = NULL,
  indices = NULL,
  n_train = 20L,
  n_features = 10L,
  n_test = 0L,
  seed = NULL,
  train_roi = NULL,
  test_roi = NULL
)
```

## Arguments

- train_data:

  Optional matrix of training data (observations x features).

- test_data:

  Optional matrix of test data (observations x features).

- indices:

  Optional integer feature indices. Defaults to
  `seq_len(ncol(train_data))`.

- n_train:

  Number of training observations to simulate when `train_data` is
  `NULL`.

- n_features:

  Number of features to simulate when `train_data` is `NULL`.

- n_test:

  Number of test observations to simulate when `test_data` is `NULL` and
  `n_test > 0`.

- seed:

  Optional RNG seed used only when simulating data.

- train_roi:

  Optional ROI object to include in the returned list.

- test_roi:

  Optional ROI object to include in the returned list.

## Value

A list with fields `train_data`, `test_data`, `indices`, `train_roi`,
and `test_roi`.
