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
  [`feature_rsa_rdm_vectors()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_rdm_vectors.md).
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
set.seed(1)
sample_ds <- gen_sample_dataset(c(4, 4, 4), nobs = 24, blocks = 3)
Fmat <- matrix(rnorm(24 * 6), 24, 6)
des <- feature_rsa_design(F = Fmat, labels = paste0("t", seq_len(24)),
                          max_comps = 3, block_var = sample_ds$design$block_var)
mdl <- feature_rsa_model(sample_ds$dataset, des, method = "pca",
                         ncomp_selection = "max", return_rdm_vectors = TRUE)
region_mask <- neuroim2::NeuroVol(
  rep(1:3, length.out = length(sample_ds$dataset$mask)),
  neuroim2::space(sample_ds$dataset$mask)
)
res <- run_regional(mdl, region_mask)
#> INFO [2026-10-08 02:06:13] 
#> MVPA Iteration Complete
#> - Total ROIs: 3
#> - Processed: 3
#> - Skipped: 0
#> INFO [2026-10-08 02:06:13] run_regional: 3 ROIs processed (success=3, errors=0)
cross_conn <- feature_rsa_cross_connectivity(res, method = "spearman")
cross_dc <- feature_rsa_cross_connectivity(
  res,
  method = "spearman",
  adjust = "double_center"
)
round(cross_dc, 2)
#>          observed
#> predicted     1     2     3
#>         1 -0.05 -0.02  0.07
#>         2 -0.02  0.04 -0.03
#>         3  0.07 -0.02 -0.04
```
