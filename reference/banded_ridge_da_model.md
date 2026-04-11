# Grouped (banded) ridge domain-adaptation model (continuous predictors -\> brain)

Fit multivariate linear models from continuous predictors (TR x
features) to ROI/searchlight activity (TR x voxels), with predictors
organized into named *feature sets*. This supports both:

- **domain adaptation**: allow or enforce similarity between
  train/source and test/target mappings, and

- **feature-set competition**: quantify what each feature set adds
  beyond the others.

Preferred name for \`banded_ridge_da_model()\`. See that function for
full details.

## Usage

``` r
banded_ridge_da_model(
  dataset,
  design,
  mode = c("stacked", "coupled"),
  lambdas,
  alpha_recall = 0.2,
  alpha_target = NULL,
  rho = 5,
  recall_folds = NULL,
  target_folds = NULL,
  recall_nfolds = 5L,
  target_nfolds = NULL,
  recall_gap = 0L,
  target_gap = NULL,
  target_nperm = 0L,
  target_perm_strategy = c("circular_shift", "block_shuffle"),
  target_perm_block = NULL,
  compute_delta_r2 = TRUE,
  delta_sets = NULL,
  return_diagnostics = FALSE,
  ...
)

grouped_ridge_da_model(
  dataset,
  design,
  mode = c("stacked", "coupled"),
  lambdas,
  alpha_recall = 0.2,
  alpha_target = NULL,
  rho = 5,
  recall_folds = NULL,
  target_folds = NULL,
  recall_nfolds = 5L,
  target_nfolds = NULL,
  recall_gap = 0L,
  target_gap = NULL,
  target_nperm = 0L,
  target_perm_strategy = "circular_shift",
  target_perm_block = NULL,
  compute_delta_r2 = TRUE,
  delta_sets = NULL,
  return_diagnostics = FALSE,
  ...
)
```

## Arguments

- dataset:

  An \`mvpa_dataset\` with \`train_data\` (train/source) and
  \`test_data\` (test/target).

- design:

  A \`feature_sets_design\` created by \`feature_sets_design()\`,
  carrying either fixed \`X_test\` predictors or a fold-aware
  \`target_builder\`.

- mode:

  \`"stacked"\` or \`"coupled"\`.

- lambdas:

  Named numeric vector of ridge penalties per feature set. Names must
  match \`names(design\$X_train\$indices)\`.

- alpha_recall:

  Non-negative scalar weighting the test/target loss (default 0.2). This
  argument is kept for backward compatibility with the original
  encoding/recall use case.

- alpha_target:

  Optional non-negative scalar weighting the test/target loss. If
  provided, overrides \`alpha_recall\`.

- rho:

  Non-negative scalar coupling strength (only used when
  \`mode="coupled"\`).

- recall_folds:

  Optional explicit list of folds, each a list with \`train\` and
  \`test\` integer indices over test rows. If NULL, folds are derived
  from \`design\$block_var_test\` when available; otherwise, contiguous
  K-fold splits are used.

- target_folds:

  Optional explicit folds over test rows. If provided, overrides
  \`recall_folds\`.

- recall_nfolds:

  Integer number of contiguous folds to use when only one test block is
  present.

- target_nfolds:

  Integer number of contiguous folds to use when only one test block is
  present. If provided, overrides \`recall_nfolds\`.

- recall_gap:

  Non-negative integer purge gap (in TRs) used only when test/target
  data consist of a single run and contiguous folds are generated
  automatically. For each held-out test segment, \`recall_gap\` TRs on
  each side are excluded from the training set to reduce temporal
  leakage from autocorrelation/physiological noise. Default 0.

- target_gap:

  Optional non-negative integer purge gap for target folds. If provided,
  overrides \`recall_gap\`.

- target_nperm:

  Non-negative integer. If \> 0 and the target domain is a single run,
  compute a permutation null by rearranging target timing while
  preserving local temporal structure (\`target_perm_strategy\`).
  Default 0 (disabled).

- target_perm_strategy:

  Permutation strategy for single-run target nulls when \`target_nperm
  \> 0\`: \`"circular_shift"\` (preserves full run autocorrelation under
  cyclic shifts) or \`"block_shuffle"\` (preserves within-block
  structure and shuffles block order).

- target_perm_block:

  Optional positive integer block size (TRs) used only when
  \`target_perm_strategy = "block_shuffle"\`. If NULL, a heuristic block
  size is used.

- compute_delta_r2:

  Logical; if TRUE, compute leave-one-set-out unique contribution
  \\\Delta R^2\\ on held-out test (default TRUE).

- delta_sets:

  Optional character vector of set names to compute \\\Delta R^2\\ for
  (default: all sets).

- return_diagnostics:

  Logical; if TRUE, store fold-level diagnostics (fold definitions,
  per-fold R-squared/MSE, hyperparameters) in
  \`regional_mvpa_result\$fits\` when running \`run_regional()\`
  (default FALSE).

- ...:

  Additional arguments (currently unused).

## Value

A model spec of class \`banded_ridge_da_model\` compatible with
\`run_regional()\` and \`run_searchlight()\`.

A model spec of class \`banded_ridge_da_model\` compatible with
\`run_regional()\` and \`run_searchlight()\`.

## Details

**Data model.** For each ROI/searchlight, rMVPA provides:

- Train/source responses \\Y\_{train}\\ from \`dataset\$train_data\`
  (rows = TRs, columns = voxels).

- Test/target responses \\Y\_{test}\\ from \`dataset\$test_data\` (rows
  = TRs, columns = voxels).

Predictors are stored on a \`feature_sets_design\`:

- Train/source predictors \\X\_{train}\\ in \`design\$X_train\$X\`.

- Test/target predictors either as fixed \\X\_{test}\\ in
  \`design\$X_test\$X\`, or as fold-specific target predictors rebuilt
  from \`design\$target_builder\`.

**Grouped (banded) ridge.** Each predictor column belongs to a named
feature set (e.g. low/mid/high/sem). The ridge penalty is applied
per-column based on that set, i.e. a diagonal penalty \\P\\ where
entries for columns in set \\g\\ take value \\\lambda_g\\.

**Two training modes.**

- *stacked*: fit a single \\\beta\\ by stacking train and test rows
  (with test down-weighted by \`alpha_recall\`/\`alpha_target\` and
  optional test TR weights from \`design\$X_test\$row_weights\`).

- *coupled*: fit \\\beta\_{train}\\ and \\\beta\_{test}\\ with coupling
  strength \`rho\`, allowing a controlled train-\>test shift while still
  borrowing strength across domains.

**Test/target-time cross-validation (default).** If test blocks/runs are
provided (\`design\$block_var_test\`) and there are at least two unique
blocks, evaluation uses leave-one-block-out on test time. If only one
block is present, evaluation falls back to contiguous folds over test
TRs (\`recall_nfolds\`/\`target_nfolds\`). In all cases, *all train TRs
are always included in training* for each fold. This avoids transductive
evaluation when test data are used in training.

**Feature-set attribution via** \\\Delta R^2\\. When \`compute_delta_r2
= TRUE\`, the model computes leave-one-set-out unique contribution on
held-out test: \$\$\Delta R^2_g = R^2\_{\mathrm{full}} - R^2\_{-g},\$\$
where \\R^2\_{-g}\\ is obtained by refitting the model without feature
set \\g\\.

## See also

[`feature_sets`](http://bbuchsbaum.github.io/rMVPA/reference/feature_sets.md),
[`expected_features`](http://bbuchsbaum.github.io/rMVPA/reference/expected_features.md),
[`feature_sets_design`](http://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md),
[`banded_ridge_da`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md),
[`grouped_ridge_da`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Build encoding predictors and declare feature sets
fs_enc <- feature_sets(X_enc, blocks(low = 100, mid = 100, high = 100, sem = 100))

# Build recall predictors from a soft alignment posterior gamma
fs_rec <- expected_features(fs_enc, gamma, drop_null = TRUE, renormalize = FALSE)

# Create design and model spec
des <- feature_sets_design(fs_enc, fs_rec, block_var_test = recall_runs)
ms <- grouped_ridge_da_model(
  dataset = dset,
  design = des,
  mode = "coupled",
  lambdas = c(low = 10, mid = 10, high = 10, sem = 10),
  alpha_recall = 0.2,
  rho = 5,
  compute_delta_r2 = TRUE,
  return_diagnostics = TRUE
)

res <- run_regional(ms, region_mask)
} # }
```
