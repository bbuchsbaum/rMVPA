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
  alphas = 10^seq(-2, 2, length.out = 9L),
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
  memory_limit_mb = 1024
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
  candidates.

- theta_method:

  One of \`"grid"\`, \`"fixed"\`, or \`"random"\`. The default
  constructs a small, deterministic simplex grid.

- theta:

  Fixed named simplex vector/matrix for \`theta_method = "fixed"\`.

- theta_grid_points:

  Grid resolution for simplex grid candidates.

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
  predictive leave-one-band-out outer-OOF delta R2. \`NULL\` (the
  default) performs no reduced-model work.

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

  Retain the complete outer out-of-fold prediction matrix.

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

  Solver intermediate-allocation contract.

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
sample <- gen_sample_dataset(c(3, 3, 3), nobs = 24,
                             response_type = "continuous", blocks = 4)
X <- matrix(rnorm(24 * 6), 24, 6)
fs <- feature_sets(X, blocks(low = 3, semantic = 3))
design <- feature_sets_design(
  fs, block_var_train = rep(1:4, each = 6), time_series = TRUE
)
model <- banded_ridge_model(
  sample$dataset, design, outer_crossval = 4, tune_crossval = 2,
  alphas = c(0.1, 1), theta_method = "fixed",
  theta = c(low = 0.5, semantic = 0.5), target_batch_size = 9,
  delta_sets = "all"
)
result <- run_banded_ridge(model)
head(result$metrics)
#>   voxel_index response n_obs       mse correlation         r2
#> 1           1  voxel_1    24 2.2490635  0.11496185 -0.3256049
#> 2           2  voxel_2    24 0.8268850  0.27293694 -0.1366843
#> 3           3  voxel_3    24 0.7820130  0.15040717 -0.2838389
#> 4           4  voxel_4    24 0.9775494 -0.50925231 -0.5134845
#> 5           5  voxel_5    24 1.2261782  0.07083063 -0.2447747
#> 6           6  voxel_6    24 1.0972899 -0.05340696 -0.5051683
head(result$predictive_leave_one_band_out$effects)
#>   voxel_index response delta_cv_r2_low delta_cv_r2_semantic
#> 1           1  voxel_1     -0.24282444          -0.09553650
#> 2           2  voxel_2     -0.08412324           0.11898830
#> 3           3  voxel_3      0.13252334          -0.23087791
#> 4           4  voxel_4     -0.01101613          -0.31467685
#> 5           5  voxel_5     -0.14902784          -0.04875116
#> 6           6  voxel_6     -0.40296980          -0.16540710
```
