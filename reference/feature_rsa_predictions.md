# Extract Per-ROI Out-of-Fold Predictions from Feature RSA Results

Pull the out-of-fold predicted patterns (\`Yhat\`) retained by
`feature_rsa_model(..., return_predictions = TRUE)`. Each ROI keeps the
merged prediction matrix, the matching observed neural patterns, the
observation order, and the outer-fold id so later scoring cannot
accidentally treat a training target as a candidate.

## Usage

``` r
feature_rsa_predictions(x)
```

## Arguments

- x:

  A `regional_mvpa_result` returned by
  [`run_regional()`](https://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
  for a `feature_rsa_model`, or a tibble/data frame that already has
  columns `roinum` and `predicted`.

## Value

A tibble with one row per ROI and columns:

- roinum:

  ROI id.

- n_obs:

  Number of observations in the retained matrices.

- observation_index:

  List-column of the observation order used for the matrices.

- fold_id:

  List-column identifying the outer test fold for each observation.
  Identification and geometry scoring must stay inside these groups.

- voxel_index:

  List-column of the spatial indices for matrix columns.

- predicted:

  List-column of out-of-fold \`Yhat\` matrices (\`n_obs\` by
  \`n_voxels\`).

- observed:

  List-column of the matching observed pattern matrices.

## See also

[`feature_rsa_model`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md),
[`feature_rsa_rdm_vectors`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_rdm_vectors.md)

## Examples

``` r
# \donttest{
  set.seed(79)
  sample <- gen_sample_dataset(c(4, 4, 4), nobs = 24, blocks = 3)
  Fmat <- matrix(rnorm(24 * 6), 24, 6)
  des <- feature_rsa_design(
    F = Fmat,
    labels = paste0("t", seq_len(24)),
    max_comps = 3,
    block_var = sample$design$block_var
  )
  mdl <- feature_rsa_model(
    sample$dataset, des, method = "pca",
    ncomp_selection = "max",
    return_predictions = TRUE
  )
  region_mask <- neuroim2::NeuroVol(
    sample(1:2, length(sample$dataset$mask), replace = TRUE),
    neuroim2::space(sample$dataset$mask)
  )
  res <- run_regional(mdl, region_mask)
#> INFO [2026-09-03 17:21:21] 
#> MVPA Iteration Complete
#> - Total ROIs: 2
#> - Processed: 2
#> - Skipped: 0
#> INFO [2026-09-03 17:21:21] run_regional: 2 ROIs processed (success=2, errors=0)
  preds <- feature_rsa_predictions(res)
  dim(preds$predicted[[1]])
#> [1] 24 34
# }
```
