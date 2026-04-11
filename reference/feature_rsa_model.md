# Create a Feature-Based RSA Model

Creates a model for feature-based Representational Similarity Analysis
(RSA) that relates neural patterns (X) to a predefined feature space
(F).

## Usage

``` r
feature_rsa_model(
  dataset,
  design,
  method = c("pls", "pca", "glmnet"),
  crossval = NULL,
  ncomp_selection = c("loo", "max", "pve"),
  pve_threshold = 0.9,
  alpha = 0.5,
  cv_glmnet = FALSE,
  lambda = NULL,
  nperm = 0,
  permute_by = c("features", "observations"),
  save_distributions = FALSE,
  return_rdm_vectors = FALSE,
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

  :   Partial Least Squares regression predicting X from F (via
      [`pls::plsr`](https://khliland.github.io/pls/reference/mvr.html)).

  pca

  :   Principal Component Regression predicting X from PCs of F (via
      [`pls::pcr`](https://khliland.github.io/pls/reference/mvr.html)).

  glmnet

  :   Elastic net regression predicting X from F using glmnet with
      multivariate Gaussian response.

- crossval:

  Optional cross-validation specification.

- ncomp_selection:

  Character string controlling how the number of components is chosen
  for `pls` and `pca` methods. One of:

  loo

  :   (Default) Fit with leave-one-out validation and select the fewest
      components within one standard error of the minimum RMSEP
      ([`pls::selectNcomp`](https://khliland.github.io/pls/reference/selectNcomp.html),
      method `"onesigma"`).

  pve

  :   Keep the fewest components whose cumulative explained variance
      reaches `pve_threshold` of the total explained by all fitted
      components.

  max

  :   Use all `max_comps` components (legacy behaviour).

  Ignored when `method = "glmnet"`.

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

  Optional numeric value or sequence of values, only used when
  method="glmnet" and cv_glmnet=FALSE. Specifies the regularization
  parameter. If NULL (default), a sequence will be automatically
  determined by glmnet.

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
  vector in the regional result's \`fits\` slot. This is off by default
  because it can add substantial memory use for long time series or many
  ROIs.

- ...:

  Additional arguments (currently unused). Passing deprecated arguments
  such as `cache_pca` now results in an error.

## Value

A `feature_rsa_model` object (S3 class).

## Details

Feature RSA models analyze how well a feature matrix `F` (defined in the
\`design\`) relates to neural data `X`. The \`max_comps\` parameter,
inherited from the \`design\` object, sets an upper limit on the number
of components fitted: - **pls**: PLS regression via
[`pls::plsr`](https://khliland.github.io/pls/reference/mvr.html). Fits
up to \`max_comps\` components; the actual number used for prediction is
chosen by `ncomp_selection`. - **pca**: Principal Component Regression
via [`pls::pcr`](https://khliland.github.io/pls/reference/mvr.html).
Fits up to \`max_comps\` components; selection controlled by
`ncomp_selection`. - **glmnet**: Elastic net regression via `glmnet`
with multivariate Gaussian response. Regularisation (lambda) can be
auto-selected via `cv_glmnet=TRUE`.

For `pls` and `pca`, the `ncomp_selection` argument determines how many
of the fitted components are actually used for prediction. The default
(`"loo"`) fits the model with leave-one-out cross-validation and picks
the fewest components within one SE of the minimum RMSEP.

\*\*Performance Metrics\*\* (computed by \`evaluate_model\` after
cross-validation):

\*Condition-pattern metrics\* (trial x trial correlation matrix): -
\`pattern_correlation\`: Average correlation between the predicted and
observed spatial patterns for corresponding trials (diagonal of the
trial x trial correlation matrix computed across voxels). -
\`pattern_discrimination\`: \`pattern_correlation\` minus the mean
off-diagonal correlation. Measures how much better the model predicts
the correct trial's pattern compared to incorrect trials. -
\`pattern_rank_percentile\`: For each trial, percentile rank of the
correct pattern match among all candidates. 0.5 = chance, 1 = perfect.

\*Representational geometry\*: - \`rdm_correlation\`: Spearman
correlation between the upper triangles of the observed and predicted
RDMs (defined as 1 - trial-by-trial correlation across voxels). Captures
similarity of representational geometry.

\*Global reconstruction metrics\*: - \`voxel_correlation\`: Correlation
of the flattened predicted and observed matrices (all trials x all
voxels). - \`mse\`: Mean Squared Error. - \`r_squared\`: 1 - RSS/TSS.

\*Voxel encoding fidelity\*: - \`mean_voxelwise_temporal_cor\`: Average
per-voxel temporal correlation between predicted and observed time
courses.

\- \`p\_\*\`, \`z\_\*\`: If \`nperm \> 0\`, permutation-based p-values
and z-scores for the above metrics.

The number of components actually used (\`ncomp\`) for the
region/searchlight is also included in the performance output.

## Examples

``` r
# \donttest{
  S <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
  labels <- factor(letters[1:5])
  des <- feature_rsa_design(S = S, labels = labels)
  # mdl <- feature_rsa_model(dataset, des, method="pls")
# }
```
