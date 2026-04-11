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
if (FALSE) { # \dontrun{
  ds <- gen_sample_dataset(c(5,5,5), 20, external_test=TRUE)
  model <- subspace_alignment_model(ds$dataset, ds$design, d=10)
} # }
```
