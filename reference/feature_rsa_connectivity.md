# Compute ROI-by-ROI Representational Connectivity from Feature RSA Predictions

Forms an ROI x ROI similarity matrix by correlating lower-triangle
predicted RDM vectors across ROIs. Sparsification, when requested, is
applied only to the final ROI x ROI matrix and never to the per-ROI RDM
vectors themselves.

## Usage

``` r
feature_rsa_connectivity(
  x,
  method = c("spearman", "pearson"),
  keep = 1,
  absolute = FALSE,
  use = "pairwise.complete.obs",
  verbose = FALSE
)
```

## Arguments

- x:

  Either a `regional_mvpa_result` produced by
  `feature_rsa_model(..., return_rdm_vectors=TRUE)` or the tibble
  returned by
  [`feature_rsa_rdm_vectors()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_rdm_vectors.md).
  Regional results may store RDM vectors either in-memory or in
  file-backed batches written by
  `run_regional(..., save_rdm_vectors_dir = ...)`.

- method:

  Correlation method used across ROI RDM vectors, one of `"spearman"` or
  `"pearson"`.

- keep:

  Proportion of ROI-ROI edges to retain after optional sparsification.
  `keep = 1` disables sparsification. For example, `keep = 0.1` retains
  the top 10% of finite off-diagonal edges.

- absolute:

  Logical; when `TRUE`, rank edges by absolute magnitude during
  sparsification. Defaults to `FALSE`.

- use:

  Missing-value handling passed to
  [`cor`](https://rdrr.io/r/stats/cor.html).

- verbose:

  Logical; if `TRUE`, emit block-level progress messages while
  connectivity is being computed.

## Value

A symmetric numeric matrix with ROIs in rows/columns.

## Examples

``` r
# Usually `vecs` comes from feature_rsa_rdm_vectors(run_regional(...));
# here a toy table of per-ROI lower-triangle RDM vectors is used.
set.seed(1)
vecs <- tibble::tibble(
  roinum = 1:4,
  rdm_vec = replicate(4, as.numeric(dist(matrix(rnorm(30), 10, 3))),
                      simplify = FALSE)
)
conn <- feature_rsa_connectivity(vecs, method = "spearman", keep = 0.5)
round(conn, 2)
#>      1     2     3     4
#> 1 1.00  0.00  0.14  0.00
#> 2 0.00  1.00 -0.02 -0.02
#> 3 0.14 -0.02  1.00  0.00
#> 4 0.00 -0.02  0.00  1.00
```
