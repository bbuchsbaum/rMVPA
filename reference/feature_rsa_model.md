# Create a Feature-Based RSA Model

Creates a model for feature-based Representational Similarity Analysis
(RSA) that relates neural patterns (X) to a predefined feature space
(F).

## Usage

``` r
feature_rsa_model(
  dataset,
  design,
  method = c("pls", "pca", "glmnet", "ridge"),
  crossval = NULL,
  ncomp_selection = c("loo", "max", "pve", "blocked"),
  pve_threshold = 0.9,
  alpha = 0.5,
  cv_glmnet = FALSE,
  lambda = NULL,
  nperm = 0,
  permute_by = c("features", "observations"),
  save_distributions = FALSE,
  return_rdm_vectors = FALSE,
  lambda_selection = c("gcv", "loo", "blocked", "fixed"),
  lambda_objective = c("mse", "pattern_discrimination", "pattern_rank_percentile"),
  lambda_one_se = NULL,
  ncomp_objective = c("mse", "pattern_discrimination", "pattern_rank_percentile"),
  ncomp_one_se = NULL,
  return_predictions = FALSE,
  max_retained_mb = 1024,
  prediction_overflow = c("error", "none"),
  feature_standardize = NULL,
  ...
)
```

## Arguments

- dataset:

  An `mvpa_dataset` object containing the neural data (`X`).

- design:

  A `feature_rsa_design` object specifying the feature space (`F`) and
  including the component limit (\`max_comps\`).

- method:

  Character string specifying the analysis method. One of:

  pls

  :   Partial Least Squares regression predicting X from F using the
      numerical algorithm configured by
      [`pls::pls.options()`](https://khliland.github.io/pls/reference/pls.options.html).

  pca

  :   Principal Component Regression predicting X from PCs of F using
      the algorithm configured by
      [`pls::pls.options()`](https://khliland.github.io/pls/reference/pls.options.html)
      (SVD-PCR by default).

  ridge

  :   Multi-response ridge regression predicting X from F with an
      economy-SVD solver and a normalized penalty.

  glmnet

  :   Elastic net regression predicting X from F using glmnet with
      multivariate Gaussian response.

- crossval:

  Optional cross-validation specification.

- ncomp_selection:

  Character string controlling how the number of components is chosen
  for `pls` and `pca` methods. One of:

  loo

  :   (Default) Use leave-one-observation-out validation and select the
      fewest components within one standard error of the minimum
      segment-wise MSE. This preserves the historical selection rule
      while streaming held-out errors instead of retaining a validation
      cube.

  pve

  :   Keep the fewest components whose cumulative explained variance
      reaches `pve_threshold` of the total explained by all fitted
      components.

  blocked

  :   Use leave-one-block-out validation within each outer training fold
      and apply the same one-standard-error rule as `"loo"`. Each
      held-out block contributes one segment MSE. Requires
      `design$block_var` and at least two training blocks in every outer
      fold.

  max

  :   Use all `max_comps` components (legacy behaviour).

  Ignored when `method` is `"ridge"` or `"glmnet"`.

- pve_threshold:

  Numeric in (0, 1\]. When `ncomp_selection = "pve"`, the proportion of
  total explained X-variance at which to stop adding components. Default
  0.9.

- alpha:

  Numeric value between 0 and 1, only used when method="glmnet".
  Controls the elastic net mixing parameter: 1 for lasso (default), 0
  for ridge, values in between for a mixture. Defaults to 0.5 (equal mix
  of ridge and lasso).

- cv_glmnet:

  Logical, if TRUE and method="glmnet", use cv.glmnet to automatically
  select the optimal lambda value via cross-validation. Defaults to
  FALSE.

- lambda:

  Optional numeric value or sequence of regularization values. For
  `method = "ridge"`, values must be finite and non-negative and
  correspond to the normalized objective \\n^{-1} \|\|Y-XB\|\|\_F^2 +
  \lambda \|\|B\|\|\_F^2\\; `NULL` uses a fixed, data-independent grid.
  For `method = "glmnet"` with `cv_glmnet = FALSE`, values must be
  positive; `NULL` lets glmnet construct its path.

- nperm:

  Integer, number of permutations to run for statistical testing of
  model performance metrics after merging cross-validation folds.
  Default 0 (no permutation testing).

- permute_by:

  DEPRECATED. Permutation is always done by shuffling rows of the
  predicted matrix.

- save_distributions:

  Logical, if TRUE and nperm \> 0, save the full null distributions from
  the permutation test. Defaults to FALSE.

- return_rdm_vectors:

  Logical; if TRUE, retain each ROI's predicted lower-triangle RDM
  vector in the regional result's \`fits\` slot. Cross-fold pairs are
  stored as missing because the two observations were not withheld
  together; the retained diagnostics include the observation order and
  fold assignment. This is off by default because it can add substantial
  memory use for long time series or many ROIs.

- lambda_selection:

  Character string controlling ridge penalty selection. `"gcv"`
  (default) minimizes generalized cross-validation error from one
  full-fold SVD. `"loo"` uses exact analytic leave-one-observation-out
  PRESS errors. `"blocked"` uses leave-one-block-out errors, with
  centering and scaling re-estimated inside every inner split; it
  requires `design$block_var`. `"fixed"` requires one supplied `lambda`.
  Selection applies the rule configured by `lambda_one_se`. Ignored for
  other methods.

- lambda_objective:

  Character string naming the ridge tuning estimand. `"mse"` (default)
  preserves GCV/LOO/blocked prediction-error tuning.
  `"pattern_discrimination"` maximizes the held-out correct-pattern
  correlation minus the mean correlation with incorrect patterns.
  `"pattern_rank_percentile"` maximizes identification rank among
  observations withheld together. Both pattern-relative objectives
  require `lambda_selection = "blocked"` with at least two observations
  per inner validation block. Ignored for other methods.

- lambda_one_se:

  Optional logical controlling the one-standard-error rule for ridge
  tuning. `NULL` uses `TRUE` for MSE-based LOO or blocked tuning and
  `FALSE` otherwise. Pattern-relative tuning defaults to the empirical
  optimum because stronger shrinkage is not inherently the safer error
  for a relative-pattern objective.

- ncomp_objective:

  Character string naming the component-count tuning estimand for `pls`
  and `pca`. `"mse"` remains the default for backward compatibility.
  `"pattern_discrimination"` maximizes the held-out correct-pattern
  correlation minus the mean correlation with incorrect patterns.
  `"pattern_rank_percentile"` maximizes held-out identification rank.
  Both pattern-relative objectives require `ncomp_selection = "blocked"`
  with at least two observations per inner validation block. Ignored for
  other methods.

- ncomp_one_se:

  Optional logical controlling the one-standard-error rule for
  component-count tuning. `NULL` uses `TRUE` for MSE-based LOO or
  blocked tuning and `FALSE` otherwise. Pattern-relative tuning defaults
  to the empirical optimum because a simpler component model is not
  automatically preferable within one standard error on these
  objectives. Ignored for other methods.

- return_predictions:

  Logical; if TRUE, retain each ROI's out-of-fold predicted patterns
  (\`Yhat\`) together with the matching observed patterns, observation
  order, outer-fold id, and voxel indices. These are stored in the
  regional result's \`fits\` slot and extracted with
  [`feature_rsa_predictions`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_predictions.md).
  Off by default because \`n_obs x n_voxels\` matrices across many ROIs
  can be large.

- max_retained_mb:

  Allocation contract (MiB) for retained out-of-fold predicted and
  observed patterns. The estimate counts both matrices for a partition
  of the active mask (one copy of each voxel). A refusal is preferred to
  a silent out-of-memory failure.

- prediction_overflow:

  What to do when the estimate exceeds `max_retained_mb`: `"error"`
  refuses the request; `"none"` disables prediction retention and
  records a notice.

- feature_standardize:

  How the feature matrix `F` is standardized inside each training fold:
  `"scale"` centers every column and divides it by its training-fold
  standard deviation; `"center"` only subtracts the training-fold column
  means. `NULL` (the default) resolves to `"center"` for designs built
  from a similarity matrix `S`, whose feature matrix consists of
  PCA-like scores, and to `"scale"` otherwise. Use `"center"` whenever
  the column variances of `F` are meaningful, e.g. PCA scores or other
  pre-whitened inputs. See the Feature standardization section.

- ...:

  Additional arguments (currently unused). Passing deprecated arguments
  such as `cache_pca` now results in an error.

## Value

A `feature_rsa_model` object (S3 class).

## Details

Feature RSA models analyze how well a feature matrix `F` (defined in the
\`design\`) relates to neural data `X`. The \`max_comps\` parameter,
inherited from the \`design\` object, sets an upper limit on the number
of components fitted: - **pls**: PLS regression using the configured pls
numerical kernel. Fits up to \`max_comps\` components, capped at the
numerical rank of the training-fold feature matrix; the actual number
used for prediction is chosen by `ncomp_selection`. - **pca**: Principal
Component Regression using the configured pls PCR kernel (SVD-PCR by
default). Fits up to \`max_comps\` components, capped in the same way;
selection is controlled by `ncomp_selection`. - **ridge**:
Multi-response ridge regression using one economy SVD per training
matrix. This is a distinct estimator, not an approximation to PLS or
elastic net. Its penalty is selected by `lambda_selection`. -
**glmnet**: Elastic net regression via `glmnet` with multivariate
Gaussian response. Regularisation (lambda) can be auto-selected via
`cv_glmnet=TRUE`.

For `pls` and `pca`, the `ncomp_selection` argument determines how many
of the fitted components are actually used for prediction. The default
(`"loo"`) uses leave-one-observation-out validation and picks the fewest
components within one SE of the minimum segment-wise MSE. The
`ncomp_objective = "mse"` default is retained for backward
compatibility. With `ncomp_selection = "blocked"`, callers may instead
maximize `"pattern_discrimination"` or `"pattern_rank_percentile"`;
these relative objectives default to the empirical optimum rather than
the one-SE rule. `"blocked"` uses leave-one-block-out segments within
each outer training fold, and centering and scaling are estimated again
from each inner training split. It is usually the more faithful and much
less expensive validation unit when observations are dependent within
acquisition runs, sessions, or subjects. It is not interchangeable with
`"pve"` or `"max"`, which do not estimate held-out performance for
component selection.

\*\*Performance Metrics\*\* (computed by \`evaluate_model\` after
cross-validation):

\*Condition-pattern metrics\* (trial x trial correlation matrix): -
\`pattern_correlation\`: Average correlation between the predicted and
observed spatial patterns for corresponding trials (diagonal of the
trial x trial correlation matrix computed across voxels). -
\`pattern_discrimination\`: \`pattern_correlation\` minus the mean
off-diagonal correlation among candidates withheld in the same outer
fold. Measures how much better the model predicts the correct trial's
pattern than eligible incorrect trials. - \`pattern_rank_percentile\`:
For each trial, percentile rank of the correct pattern match among
candidates withheld in the same outer fold. 0.5 = chance, 1 = perfect.

\*Representational geometry\*: - \`rdm_correlation\`: Spearman
correlation between jointly withheld pairs in the observed and predicted
RDMs (defined as 1 - trial-by-trial correlation across voxels). Captures
similarity of held-out representational geometry without comparing a
prediction with a target that trained that prediction.

\*Global reconstruction metrics\*: - \`voxel_correlation\`: Correlation
of the flattened predicted and observed matrices (all trials x all
voxels). - \`mse\`: Mean Squared Error. - \`r_squared\`: 1 - RSS/TSS.

\*Voxel encoding fidelity\*: - \`mean_voxelwise_temporal_cor\`: Average
per-voxel temporal correlation between predicted and observed time
courses.

\- \`p\_\*\`, \`z\_\*\`: If \`nperm \> 0\`, permutation-based p-values
and z-scores for the above metrics.

For PLS, PCR, and glmnet, the number of components or its historical
proxy (\`ncomp\`) is included in the performance output. Ridge instead
reports the median selected penalty (\`median_lambda\`), mean effective
degrees of freedom (\`mean_effective_df\`), non-intercept degrees of
freedom, and fractions of folds selected at either end of the lambda
grid.

\*\*Out-of-fold predictions\*\* (\`return_predictions = TRUE\`): Each
ROI retains the merged out-of-fold \`Yhat\` and \`Y\` matrices, the
observation order, and the outer-fold id used by the built-in scoring
rules. Extract them with
[`feature_rsa_predictions`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_predictions.md)
for post-hoc scoring (temporal windowing, whitened distances, custom
identification rules) without crossing fold boundaries. This is a
different payload from \`return_rdm_vectors\` and from the
classification \`prediction_table\`. Overlapping searchlights are
refused because they would retain multiple copies of each voxel.

## Feature standardization

Inside every training fold the neural responses `X` are always z-scored
column-wise, and by default the feature matrix `F` is too
(`feature_standardize = "scale"`). Held-out rows are transformed with
the training-fold means and scales, and predictions are returned on the
original response scale. This differs from the default of
[`pls::pcr()`](https://khliland.github.io/pls/reference/mvr.html)
(`scale = NULL`), so `method = "pca"` is correlation-PCR rather than
covariance-PCR.

Column scaling discards the variance profile of `F`. That is harmless
for raw feature spaces but degenerate for inputs whose column variances
carry the structure, such as PCA scores or any whitened matrix: after
scaling every direction has unit variance, PCA component ordering is set
by noise, and `F = U %*% S` gives exactly the same fit as `F = U`. Use
`feature_standardize = "center"` for such inputs; it centers `F` per
training fold without rescaling columns, and applies to the PLS/PCA,
ridge, and glmnet paths alike, including nested blocked tuning. Designs
built from a similarity matrix `S` store exactly such scores
(eigenvectors weighted by the square roots of the eigenvalues), so they
default to `"center"`.

For `method = "pca"` under column scaling, the constructor inspects the
correlation spectrum of `F` and warns when it is flat (condition number
within two percent of one), the signature of columns that are exactly
orthogonal on these rows, such as PCA scores computed from them. Scores
computed on a superset of rows, or otherwise nearly whitened inputs, are
not detected, so use `"center"` for any PCA-score input whether or not a
warning appears. The diagnostic is stored in `model$feature_spectrum`.

Degenerate columns of `F` are treated according to what the requested
standardization actually does. Under `"scale"` a column that is
constant, or near-constant relative to the widest column, cannot be
divided by its standard deviation: a column that is constant over every
row is refused here, by the constructor, and one that is constant only
within some training fold is refused by that fit, in both cases with a
message naming how many such columns there are and which is worst.
Constancy is judged against each column's own magnitude, never against
the other columns' or against an absolute threshold, so `F` may mix
units freely and rescaling it never changes the outcome. Under
`"center"` no division takes place, so those columns are accepted; for
`method = "pls"` and `method = "pca"` the number of components is then
capped at the numerical rank of the feature matrix, both for the final
fit and within each segment of blocked or leave-one-out component
selection, which is what keeps undefined directions out of the fit and
out of the tuning scores.

## See also

[`feature_rsa_predictions`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_predictions.md),
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
#> INFO [2026-09-04 17:44:14] 
#> MVPA Iteration Complete
#> - Total ROIs: 2
#> - Processed: 2
#> - Skipped: 0
#> INFO [2026-09-04 17:44:15] run_regional: 2 ROIs processed (success=2, errors=0)
  preds <- feature_rsa_predictions(res)
  dim(preds$predicted[[1]])
#> [1] 24 34
# }
```
