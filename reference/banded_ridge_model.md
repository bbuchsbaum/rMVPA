# First-class single-domain banded-ridge encoding model

\`banded_ridge_model()\` defines leakage-safe nested tuning for grouped
stimulus predictors and independent brain responses. Feature-band
membership and training blocks come exclusively from
\`feature_sets_design()\`.

## Usage

``` r
banded_ridge_model(
  dataset,
  design,
  outer_crossval = NULL,
  tune_crossval = NULL,
  candidates = NULL,
  alphas = "auto",
  theta_method = c("grid", "fixed", "random"),
  theta = NULL,
  theta_grid_points = 3L,
  n_theta = 20L,
  metric = c("mse", "correlation", "r2"),
  alpha_scope = "response",
  theta_scope = "response",
  fixed_alpha = NULL,
  fixed_theta = NULL,
  delta_sets = NULL,
  purge = 0L,
  solver = c("auto", "direct", "svd_primal", "dual_kernel"),
  target_batch_size = 256L,
  return_predictions = FALSE,
  weight_retention = c("none", "dual", "primal"),
  retain_diagnostics = FALSE,
  max_retained_mb = 1024,
  weight_overflow = c("error", "none"),
  seed = 1L,
  memory_limit_mb = 1024,
  nuisance = NULL,
  nuisance_fit = c("train", "within", "global")
)
```

## Arguments

- dataset:

  An image \`mvpa_dataset\` whose training rows are observations and
  active mask locations are responses.

- design:

  A \`feature_sets_design\` containing \`X_train\`, band membership, and
  optional \`block_var_train\` / \`time_series\` metadata.

- outer_crossval:

  Outer fold specification: an integer fold count, a
  \`cross_validation\` object (for example
  [`blocked_cross_validation`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md),
  [`kfold_cross_validation`](https://bbuchsbaum.github.io/rMVPA/reference/kfold_cross_validation.md),
  or
  [`custom_cross_validation`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md);
  only schemes whose test partitions cover every training row exactly
  once convert, so resampled or non-deterministic schemes are refused),
  or an explicit list of \`list(id =, train =, test =)\` folds. If
  omitted, intact training blocks are required and define at most five
  folds. Integer counts and \`cross_validation\` objects have \`purge\`
  applied to their training rows; an explicit fold list is validated as
  supplied and must already respect \`purge\`.

- tune_crossval:

  Inner fold count, or one explicit inner-fold list per outer fold.
  \`NULL\` (the default) uses at most five inner folds, limited by the
  number of declared training blocks (or rows, when no blocks are
  declared) available inside each outer training set. A supplied value
  is still type-checked but otherwise ignored when only one candidate
  survives the alpha/theta scopes (for example \`theta_method =
  "fixed"\` with one \`theta\` and one \`alphas\` value, or fixed scopes
  with \`fixed_alpha\` and \`fixed_theta\`): nothing needs selecting, so
  inner tuning is skipped, the model's \`inner_v\` is \`NA\`, and
  reported \`inner_score\` values are \`NA\`. Inner folds are
  constructed and validated when the model is created, so an impossible
  request fails here rather than in \`run_banded_ridge()\`.

- candidates:

  Optional candidate manifest. By default one is created from
  \`alphas\`, \`theta_method\`, and the remaining theta arguments.

- alphas:

  Non-negative overall ridge penalties used when constructing
  candidates, or \`"auto"\` (the default) to scale the grid to the
  design. Every solver path standardizes each training column and then
  scales band \`b\` by \`sqrt(theta_b)\`, so \`alpha\` competes with the
  eigenvalues of the resulting cross-product, which average \`(n - 1) \*
  mean(D_b) / min(n - 1, p)\` under a uniform \`theta\` for band widths
  \`D_b\`. \`"auto"\` places nine points spanning six decades centred on
  that anchor, so the grid brackets the range in which the penalty
  changes the fit however many rows and features the design has. A fixed
  grid cannot: a design with a few thousand timepoints needs penalties
  several orders of magnitude above the fixed \`10^seq(-2, 2)\` grid
  this default replaced, and every response then selects the largest
  available value. \`run_banded_ridge()\` warns when that happens; see
  \`result\$selection_diagnostics\`. The resolved grid is on the
  returned model as \`alpha_grid\`.

- theta_method:

  One of \`"grid"\`, \`"fixed"\`, or \`"random"\`. The default
  constructs a small, deterministic simplex grid.

- theta:

  Fixed named simplex vector/matrix for \`theta_method = "fixed"\`.

- theta_grid_points:

  Grid resolution for simplex grid candidates. A grid with fewer than
  \`n_bands + 1\` points has no candidate in which every band is active,
  so the equal-weight point is added automatically. With more than about
  three bands prefer \`theta_method = "random"\`, which samples the
  interior of the simplex.

- n_theta:

  Number of deterministic random simplex points.

- metric:

  Inner selection metric; MSE is the default estimand.

- alpha_scope, theta_scope:

  \`"response"\`, \`"shared"\`/\`"roi"\`, or \`"fixed"\` selection
  constraints.

- fixed_alpha, fixed_theta:

  Values required by corresponding fixed scopes.

- delta_sets:

  Optional feature-band names, or \`"all"\`, for independently retuned
  predictive leave-one-band-out outer-OOF delta R2. A named list of
  band-name vectors drops each entry's bands together, for example
  \`list(phonology = c("manner", "place", "voicing"), manner =
  "manner")\`, giving the unique contribution of a family of bands next
  to that of its members; results are labelled by the list names.
  \`"families"\` requests every band plus every family declared in
  \`feature_sets(families = ...)\`. \`NULL\` (the default) performs no
  reduced-model work.

- purge:

  Non-negative temporal gap (in rows) around validation/test rows.
  Training rows within the gap are dropped when folds are constructed
  from an integer count or a \`cross_validation\` object; explicit fold
  lists (outer or inner) are checked against the gap and rejected if
  they violate it, never modified.

- solver:

  \`"auto"\`, \`"direct"\`, \`"svd_primal"\`, or \`"dual_kernel"\`.

- target_batch_size:

  Maximum responses evaluated together.

- return_predictions:

  Retain the complete outer out-of-fold prediction matrix. With nuisance
  regressors, the nuisance-residualized held-out responses the
  predictions are scored against are also returned, as \`observed\`.

- weight_retention:

  Retain no weights, outer-fold primal coefficients, or outer-fold dual
  weights.

- retain_diagnostics:

  Retain full candidate/fold diagnostics per chunk.

- max_retained_mb:

  Allocation contract for predictions, weights, and diagnostics.

- weight_overflow:

  Refuse an unsafe request or fall back to no weights.

- seed:

  Deterministic fold/candidate seed.

- memory_limit_mb:

  Solver intermediate-allocation contract in MiB. It bounds two separate
  pools, so peak solver memory can approach twice this value: per-fit
  intermediates are checked against it up front, and the retained
  decomposition cache of the optimized solvers is capped at it
  independently. Once cached band Grams and eigendecompositions would
  exceed the cap, the least recently used entries are evicted and
  recomputed on demand. Results are unchanged at any value; only reuse
  is lost, and provenance reports \`solver_cache_evictions\` and
  \`solver_cache_peak_mb\`.

- nuisance:

  Optional numeric matrix of nuisance regressors (rows = observations),
  for example motion parameters and drift terms. Within every inner and
  outer split an intercept and these columns are regressed out of both
  the features and the responses, separately for the training rows and
  the held-out rows (see \`nuisance_fit\`). Nuisance columns are never
  penalized and never contribute to the predictions, so all reported
  accuracies, including leave-band-out delta R2, describe variance in
  the nuisance-residualized held-out responses.

- nuisance_fit:

  How held-out rows are residualized. \`"train"\` (default) applies the
  coefficients fitted on the training rows, so the transformation is
  fixed before the held-out rows are seen. This requires the held-out
  nuisance rows to lie within the range of the training rows. With
  drift-like regressors and contiguous held-out blocks they do not, the
  extrapolated fit would dominate features and responses alike and
  inflate accuracy, and the model stops with an error. \`"within"\` fits
  the nuisance coefficients on the held-out rows themselves. Nothing is
  carried between training and held-out data and nothing is
  extrapolated, but a held-out row's transformation then depends on the
  other rows of its assessment set, so the reported accuracy is
  conditional on each held-out block and is not the accuracy of a frozen
  training-derived pipeline. Each held-out split needs more rows than
  nuisance regressors. \`"global"\` removes the nuisance from features
  and responses once, over all rows, before the data are split, as much
  encoding work does in preprocessing. \*\*It is not leak-free and its
  accuracies can be biased upward.\*\* Held-out responses enter the
  nuisance coefficients that are then subtracted from the training
  responses, so training and held-out residuals are dependent, and
  responses with no relation to the features can show positive held-out
  accuracy. The bias grows with the number of nuisance and feature
  columns relative to the number of rows: in a 374-volume run with 25
  nuisance and 44 feature columns, unrelated responses gave held-out
  correlations of about .03 to .05 where \`"within"\` gave about zero,
  and when few residual degrees of freedom remain the fit becomes exact.
  Use it for exploration, or calibrate it against a null computed the
  same way (for example with time-shifted features); do not test its
  accuracies against zero. The model refuses when the features have as
  many columns as there are residual degrees of freedom. The mode used
  is recorded in \`result\$provenance\$nuisance_fit\`.

## Value

A \`banded_ridge_model\` specification for \`run_banded_ridge()\`.

## Predictive leave-one-band-out delta R2

When \`delta_sets\` is enabled, each requested reduced model is tuned
again inside every outer-training set, over a candidate manifest
projected onto the retained bands. The reported effect is \`R2_full -
R2_without_band\` from matched outer out-of-fold predictions. It is a
predictive drop-out effect, not an additive unique/shared variance
partition: values are not clipped, may be negative, and need not sum to
the full-model R2.

## Examples

``` r
set.seed(71)
n <- 24L
dims <- c(3L, 3L, 3L)
X <- matrix(rnorm(n * 6), n, 6)
Y <- X %*% matrix(rnorm(6 * prod(dims), sd = 0.8), 6, prod(dims)) +
  matrix(rnorm(n * prod(dims), sd = 0.5), n, prod(dims))
dataset <- mvpa_dataset(
  neuroim2::NeuroVec(array(as.vector(t(Y)), c(dims, n)),
                     neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))),
  mask = neuroim2::NeuroVol(array(1, dims),
                            neuroim2::NeuroSpace(dims, c(1, 1, 1)))
)
fs <- feature_sets(X, blocks(low = 3, semantic = 3))
design <- feature_sets_design(
  fs, block_var_train = rep(1:4, each = 6), time_series = TRUE
)
model <- banded_ridge_model(
  dataset, design, outer_crossval = 4, tune_crossval = 2,
  theta_method = "fixed", theta = c(low = 0.5, semantic = 0.5),
  target_batch_size = 9, delta_sets = "all"
)
result <- run_banded_ridge(model)
head(result$metrics)
#>   voxel_index response n_obs       mse correlation        r2
#> 1           1  voxel_1    24 0.5168105   0.9657918 0.9168180
#> 2           2  voxel_2    24 0.3340694   0.9618700 0.9229610
#> 3           3  voxel_3    24 0.7278631   0.7454001 0.5520757
#> 4           4  voxel_4    24 0.5893974   0.9316744 0.8622335
#> 5           5  voxel_5    24 0.2821740   0.7864845 0.6128238
#> 6           6  voxel_6    24 0.3863853   0.9233786 0.8493119
result$selection_diagnostics$alpha
#>              model n_selections n_alpha_grid alpha_grid_min alpha_grid_max
#> 1             full          108            9         0.0115          11500
#> 2      without_low          108            9         0.0115          11500
#> 3 without_semantic          108            9         0.0115          11500
#>   modal_alpha modal_share share_at_grid_min share_at_grid_max share_interior
#> 1   0.3636619   0.3333333        0.24074074        0.00000000      0.7592593
#> 2   2.0450213   0.3518519        0.07407407        0.14814815      0.7777778
#> 3   2.0450213   0.2407407        0.20370370        0.06481481      0.7314815
#>   saturated
#> 1     FALSE
#> 2     FALSE
#> 3     FALSE
head(result$predictive_leave_one_band_out$effects)
#>   voxel_index response delta_cv_r2_low delta_cv_r2_semantic
#> 1           1  voxel_1       0.8524352           0.45597131
#> 2           2  voxel_2       0.6055967           0.19375823
#> 3           3  voxel_3       0.3164995           0.43824971
#> 4           4  voxel_4       0.2287460           0.31607592
#> 5           5  voxel_5       1.0416128           0.27988198
#> 6           6  voxel_6       1.1002373           0.06409502
```
