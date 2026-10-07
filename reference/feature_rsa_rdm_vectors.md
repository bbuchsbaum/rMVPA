# Extract Per-ROI Predicted and Observed RDM Vectors from Feature RSA Results

Convenience helper to pull the compact lower-triangle predicted (and
optionally observed) RDM vectors stored by
`feature_rsa_model(..., return_rdm_vectors = TRUE)` from a
`regional_mvpa_result`. Results can come either from in-memory `$fits`
or from file-backed batches written via
`run_regional(..., save_rdm_vectors_dir = ...)`.

## Usage

``` r
feature_rsa_rdm_vectors(x)
```

## Arguments

- x:

  A `regional_mvpa_result` returned by
  [`run_regional()`](https://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
  for a `feature_rsa_model`, or a tibble/data frame with columns
  `roinum` and `rdm_vec`. A `regional_mvpa_result` may store vectors
  either in-memory or on disk in `$rdm_batch_dir`.

## Value

A tibble with one row per ROI and columns:

- roinum:

  ROI id.

- n_obs:

  Number of observations contributing to the vector.

- observation_index:

  List-column of observation ordering used for the predicted RDM.

- fold_id:

  List-column identifying the outer test fold for each observation.
  Cross-fold entries in the RDM vectors are missing.

- rdm_vec:

  List-column containing the lower-triangle predicted RDM vector for
  that ROI.

- observed_rdm_vec:

  List-column containing the lower-triangle observed RDM vector for that
  ROI (if available).

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
#> INFO [2026-10-07 05:02:00] 
#> MVPA Iteration Complete
#> - Total ROIs: 3
#> - Processed: 3
#> - Skipped: 0
#> INFO [2026-10-07 05:02:00] run_regional: 3 ROIs processed (success=3, errors=0)
vecs <- feature_rsa_rdm_vectors(res)
vecs
#> # A tibble: 3 × 6
#>   roinum n_obs observation_index fold_id    rdm_vec     observed_rdm_vec
#>    <int> <int> <list>            <list>     <list>      <list>          
#> 1      1    24 <int [24]>        <int [24]> <dbl [276]> <dbl [276]>     
#> 2      2    24 <int [24]>        <int [24]> <dbl [276]> <dbl [276]>     
#> 3      3    24 <int [24]>        <int [24]> <dbl [276]> <dbl [276]>     
```
