# Compute Second-Order Similarity Scores

Calculates correlation-based *second order similarity* between:

- A **full NxN distance matrix** computed from `X` via `distfun`, and

- A `Dref` matrix (the "reference" dissimilarities).

For each row `i`, this excludes same-block comparisons by selecting
`which(block_var != block_var[i])`.

## Usage

``` r
second_order_similarity(
  distfun,
  X,
  Dref,
  block_var,
  method = c("pearson", "spearman")
)
```

## Arguments

- distfun:

  An S3 distance object (see
  [`create_dist`](http://bbuchsbaum.github.io/rMVPA/reference/create_dist.md))
  specifying how to compute a pairwise distance matrix from `X`.

- X:

  A numeric matrix (rows = observations, columns = features).

- Dref:

  A numeric NxN reference matrix of dissimilarities (e.g., from an ROI
  mask or a prior).

- block_var:

  A vector indicating block/group memberships for each row in `X`.

- method:

  Correlation method: "pearson" or "spearman".

## Value

A numeric vector of length `nrow(X)`, where each entry is the
correlation (using `method`) between `distance_matrix[i, valid]` and
`Dref[i, valid]`, with `valid = which(block_var != block_var[i])`.

## Details

This function first calls `pairwise_dist(distfun, X)`, obtaining an NxN
matrix of *all* pairwise distances. It does not do block-based exclusion
internally. Instead, for each row `i`, it excludes same-block rows from
the correlation by subsetting the distances to `valid_indices`.

## Examples

``` r
# Suppose we have X (10x5), a reference D (10x10), block var, and a correlation distfun:
X <- matrix(rnorm(50), 10, 5)
D <- matrix(runif(100), 10, 10)
block <- rep(1:2, each=5)
dist_obj <- cordist(method="pearson")
scores <- second_order_similarity(dist_obj, X, D, block, method="spearman")
```
