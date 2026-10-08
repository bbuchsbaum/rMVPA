# Subspace Alignment cross-decoder

Fast unsupervised domain-adaptation baseline following Fernando et al.
(ICCV 2013). Learns separate PCA subspaces for source (train) and target
(test), aligns them with a closed-form map \\M = X_S^T X_T\\, projects
both domains, and classifies target trials via correlation to source
class prototypes in the aligned space. Requires an external test set but
no target labels for fitting.

## Usage

``` r
subspace_alignment_model(
  dataset,
  design,
  d = 20L,
  center = TRUE,
  scale = TRUE,
  return_predictions = TRUE,
  ...
)
```

## Arguments

- dataset:

  mvpa_dataset with \`train_data\` (source) and \`test_data\` (target).

- design:

  mvpa_design with \`y_train\` (source labels) and \`y_test\` (for
  evaluation).

- d:

  Integer subspace dimension; capped automatically by samples/features.

- center, scale:

  Logical flags for per-domain z-normalization prior to PCA.

- return_predictions:

  logical; keep per-ROI predictions (default TRUE).

- ...:

  Additional arguments stored on the model spec.

## Value

A model spec of class \`subspace_alignment_model\` for use with
\`run_regional()\` / \`run_searchlight()\`.

## Examples

``` r
ds <- gen_sample_dataset(c(4, 4, 4), 24, nlevels = 3, blocks = 3,
                         external_test = TRUE)
model <- subspace_alignment_model(ds$dataset, ds$design, d = 5)
region_mask <- neuroim2::NeuroVol(array(1, c(4, 4, 4)),
                                  neuroim2::space(ds$dataset$mask))
res <- run_regional(model, region_mask)
#> INFO [2026-10-08 15:55:50] 
#> MVPA Iteration Complete
#> - Total ROIs: 1
#> - Processed: 1
#> - Skipped: 0
#> INFO [2026-10-08 15:55:50] run_regional: 1 ROIs processed (success=1, errors=0)
res$performance_table
#> # A tibble: 1 × 5
#>   roinum Accuracy    AUC d_used alignment_frob
#>    <int>    <dbl>  <dbl>  <dbl>          <dbl>
#> 1      1    0.292 -0.115      5           2.13
```
