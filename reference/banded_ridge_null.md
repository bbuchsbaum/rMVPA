# Permutation null for banded-ridge prediction accuracy

Tests the held-out accuracy of every response against a null in which
the predicted series is shifted in time against the observed one.
Shifting keeps the autocorrelation of both series, which a row shuffle
would destroy, so the null is appropriate for continuous time series. No
model is refit.

## Usage

``` r
banded_ridge_null(
  model,
  result,
  n_perm = 1000L,
  metric = c("correlation", "r2"),
  strategy = c("circular_shift", "block_shuffle"),
  block_size = NULL,
  min_shift = 1L,
  seed = 1L
)
```

## Arguments

- model:

  A \`banded_ridge_model\` built with \`return_predictions = TRUE\`.

- result:

  The result of \`run_banded_ridge(model)\`.

- n_perm:

  Number of permutations.

- metric:

  \`"correlation"\` (default) or \`"r2"\`.

- strategy:

  \`"circular_shift"\` (default) or \`"block_shuffle"\`.

- block_size:

  Block length for \`"block_shuffle"\`; defaults to \`floor(sqrt(n))\`.

- min_shift:

  Smallest circular shift, in rows. Set it beyond the autocorrelation
  length of the data.

- seed:

  Integer seed; the same permutations are applied to every response.

## Value

A data frame with one row per response: \`voxel_index\`, the observed
\`statistic\`, the permutation \`p_value\` (one-sided, observed value
counted in the null), and \`z\`, the statistic standardized by the null.

## Details

The null addresses the full model's accuracy. It does not test a band's
unique contribution: shifting the predictions of the full and reduced
models together leaves their difference near zero by construction.

## See also

\[run_banded_ridge()\]
