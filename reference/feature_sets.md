# Feature sets: grouped predictor matrices

Many rMVPA analyses use continuous predictors (rows = TRs/observations,
columns = stimulus features) to explain continuous neural responses
(rows = TRs, columns = voxels/parcels). In practice, stimulus predictors
often come in *multiple correlated blocks* (e.g. VGG
low/mid/high/semantic PCs, or audio vs vision features), and it is
useful to:

- regularize each block differently (banded/grouped ridge), and

- compute attribution/competition measures such as leave-one-set-out
  \\\Delta R^2\\.

Wrap a predictor representation into a \`feature_sets\` instance.

## Usage

``` r
feature_sets(x, spec = NULL, set_order = NULL, row_weights = NULL)
```

## Arguments

- x:

  Either (1) a numeric matrix (observations x features) with a set spec,
  or (2) a named list of numeric matrices, one per set (all with the
  same number of rows), which will be column-bound in \`set_order\`.

- spec:

  A set spec as created by \`blocks()\` or \`by_set()\` (required when
  \`x\` is a matrix).

- set_order:

  Optional character vector giving the desired set order (defaults to
  the order implied by \`spec\` or the names of the list).

- row_weights:

  Optional numeric vector of length nrow(X), used as observation weights
  by downstream models (default: all 1).

## Value

An object of class \`feature_sets\`.

## Details

The \`feature_sets()\` constructor wraps a predictor matrix (or a list
of per-set matrices) into a \`feature_sets\` object that carries set
labels and column indices so that downstream models do not need to
manually track slices.

You can supply predictors in two equivalent ways:

1.  A single matrix `X` (observations x features) plus a set
    specification via \`blocks()\` or \`by_set()\`.

2.  A named list of matrices, one per set, each with the same number of
    rows. The matrices are column-bound in \`set_order\`.

The returned object includes:

- \`X\`: the concatenated numeric predictor matrix,

- \`set\`: a factor of length \`ncol(X)\` giving per-column set
  membership,

- \`indices\`: a named list mapping set name -\> integer column indices,

- \`dims\`: number of columns per set,

- \`row_weights\`: optional observation weights (length \`nrow(X)\`).

## Feature set specifications

Use \`blocks()\` when your columns are already concatenated into
consecutive blocks (e.g. `[low|mid|high|sem]`). Use \`by_set()\` when
you have a per-column set label vector (non-contiguous groups).

## Row weights

A \`feature_sets\` object also stores \`row_weights\`. These are
optional observation weights (length = `nrow(X)`). They are primarily
intended for soft-alignment use cases where recall predictors are
computed as `gamma %*% X_enc` and the remaining posterior mass can be
used to down-weight uncertain recall TRs (see \`expected_features()\`).

## See also

[`blocks`](http://bbuchsbaum.github.io/rMVPA/reference/blocks.md),
[`by_set`](http://bbuchsbaum.github.io/rMVPA/reference/by_set.md),
`feature_sets`,
[`expected_features`](http://bbuchsbaum.github.io/rMVPA/reference/expected_features.md),
[`feature_sets_design`](http://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md),
[`banded_ridge_da_model`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md),
[`grouped_ridge_da_model`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md),
[`banded_ridge_da`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md),
[`grouped_ridge_da`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md)

[`blocks`](http://bbuchsbaum.github.io/rMVPA/reference/blocks.md),
[`by_set`](http://bbuchsbaum.github.io/rMVPA/reference/by_set.md),
[`expected_features`](http://bbuchsbaum.github.io/rMVPA/reference/expected_features.md)

## Examples

``` r
X <- matrix(rnorm(20 * 8), 20, 8)
fs <- feature_sets(X, blocks(low = 3, sem = 5))
fs
#> feature_sets
#> ===========
#> 
#> Observations: 20
#> Features:     8
#> Sets:         2 (low, sem)
# 1) Matrix input + blocks()
X <- matrix(rnorm(20 * 8), 20, 8)
fs <- feature_sets(X, blocks(low = 3, sem = 5))

# 2) List input (already split per set)
Xlist <- list(low = X[, 1:3], sem = X[, 4:8])
fs2 <- feature_sets(Xlist)
```
