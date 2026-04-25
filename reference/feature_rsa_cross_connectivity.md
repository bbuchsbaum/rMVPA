# Compute Cross-Connectivity: Predicted-Observed ROI x ROI Matrix

Builds an asymmetric ROI x ROI matrix where entry (i, j) is the
correlation between the predicted RDM vector of ROI i and the observed
RDM vector of ROI j. This captures how well the model-predicted
representational geometry in one ROI matches the data-driven geometry in
another.

## Usage

``` r
feature_rsa_cross_connectivity(
  x,
  method = c("spearman", "pearson"),
  adjust = c("none", "double_center", "residualize_mean"),
  return_components = FALSE,
  use = "pairwise.complete.obs",
  verbose = FALSE
)
```

## Arguments

- x:

  Either a `regional_mvpa_result` produced by
  `feature_rsa_model(..., return_rdm_vectors=TRUE)` or the tibble
  returned by
  [`feature_rsa_rdm_vectors()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_rdm_vectors.md).
  Regional results may store RDM vectors either in-memory or in
  file-backed batches written by
  `run_regional(..., save_rdm_vectors_dir = ...)`.

- method:

  Correlation method, one of `"spearman"` or `"pearson"`.

- adjust:

  Optional adjustment for ROI-level source/target offsets. Use `"none"`
  (default) to return the raw ROI x ROI correlation matrix,
  `"double_center"` to subtract additive source and target main effects
  from that matrix, or `"residualize_mean"` to remove the grand-mean RDM
  component from predicted and observed ROI vectors before computing the
  cross-correlation.

- return_components:

  Logical; if `TRUE`, return a list containing the requested matrix, the
  raw matrix, the adjusted matrix, and the source and target offset
  terms estimated from the raw matrix.

- use:

  Missing-value handling passed to
  [`cor`](https://rdrr.io/r/stats/cor.html).

- verbose:

  Logical; if `TRUE`, emit block-level progress messages while
  cross-connectivity is being computed.

## Value

By default, a numeric matrix of dimension n_ROI x n_ROI. Rows correspond
to predicted RDM vectors and columns to observed RDM vectors. The matrix
is *not* necessarily symmetric. If `return_components = TRUE`, a list is
returned with elements `matrix`, `raw_matrix`, `adjusted_matrix`,
`source_offset`, `target_offset`, `grand_mean`, `method`, and `adjust`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_regional(
  feature_rsa_model(dataset, design, method = "pls", return_rdm_vectors = TRUE),
  region_mask
)
cross_conn <- feature_rsa_cross_connectivity(res, method = "spearman")
cross_dc <- feature_rsa_cross_connectivity(
  res,
  method = "spearman",
  adjust = "double_center"
)
} # }
```
