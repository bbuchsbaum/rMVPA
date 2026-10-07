# Run Permutation Searchlight Inference

Computes permutation-based p-values for a searchlight MVPA result by
running the analysis on permuted labels, building a covariate-adjusted
null distribution, and mapping p-values back to all brain voxels.

## Usage

``` r
run_permutation_searchlight(
  model_spec,
  observed = NULL,
  radius = 8,
  method = c("standard", "randomized"),
  perm_ctrl = permutation_control(),
  metric = NULL,
  ...
)
```

## Arguments

- model_spec:

  An `mvpa_model`, `rsa_model`, or compatible model specification whose
  design has a
  [`permute_labels`](https://bbuchsbaum.github.io/rMVPA/reference/permute_labels.md)
  method. For `rsa_model` the permutation relabels items, so the null
  carries the same item-level dependence as the observed statistic.

- observed:

  Optional pre-computed searchlight result (output of
  [`run_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)).
  If `NULL`, it is computed internally. Can also be a named numeric
  vector of observed metric values indexed by center ID.

- radius:

  Numeric. Searchlight radius in mm (default 8).

- method:

  Character. Searchlight method passed to
  [`run_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  when `observed` is `NULL` or when `perm_strategy = "searchlight"`.

- perm_ctrl:

  A
  [`permutation_control`](https://bbuchsbaum.github.io/rMVPA/reference/permutation_control.md)
  object.

- metric:

  Character vector naming the performance metric(s) to test. If `NULL`
  (default), the first metric is used, except for `rsa_model`
  specifications, where every model predictor is tested: a permutation
  pass returns all predictors at once, so scoring them uses the same
  permutation passes. The null hypothesis is controlled by
  `perm_ctrl$rsa_null`, not by how many metrics are selected. A metric
  that is named explicitly but absent from the results is an error
  rather than a silent fallback to the first metric.

- ...:

  Additional arguments forwarded to
  [`run_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  (observed pass and `"searchlight"` permutations). For `"iterate"`
  permutations, only arguments that are formal parameters of
  [`mvpa_iterate`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
  are forwarded; other keys are ignored for that path to avoid
  argument-mismatch failures (e.g., `engine = "legacy"` is meaningful
  for `run_searchlight` but not for `mvpa_iterate`).

## Value

When one metric is tested, a `permutation_result` S3 object. When
several are tested (the default for `rsa_model`), a
`permutation_result_set`: a named list of `permutation_result` objects,
one per metric, all scored against the same permutations. Each
`permutation_result` contains:

- p_map:

  Spatial map of raw p-values.

- p_adj_map:

  Spatial map of FDR-adjusted p-values (if requested).

- p_values:

  Numeric vector of raw p-values (all centers).

- p_adjusted:

  Numeric vector of adjusted p-values.

- observed:

  The observed searchlight result.

- diagnostics:

  A `"null_diagnostics"` object (if `diagnose = TRUE`).

- perm_ctrl:

  The `permutation_control` used.

- metric:

  Metric name used for inference.

- perm_strategy:

  The strategy that was actually used.

- rsa_null, null_hypothesis, null_predictors:

  For RSA, the null scope, its description, and the predictors covered
  by that null. These fields are `NULL` for other model classes.

## RSA null hypothesis

Raw item permutations break associations with every design predictor,
including nuisance predictors. With `regtype = "lm"` or `"rfit"` and
multiple design predictors, they do not preserve the effects of the
other predictors when testing one coefficient. Such runs are refused by
default, even when `metric` selects a single output. Use
`permutation_control(rsa_null = "joint")` only to test the joint
no-association null. Each output metric then supplies a statistic for
that same joint null; significance does not establish that the
corresponding predictor has a unique effect. Conditional coefficient
inference is not implemented. This restriction also covers semi-partial
and constrained regression statistics.

Correlation models test marginal associations, without adjusting for the
other RDMs, and single-predictor regressions remain available with the
default `rsa_null = "individual"`. RSA results record `rsa_null`,
`null_hypothesis`, and `null_predictors`. FDR adjustment is performed
separately within each metric's spatial map; it does not correct across
metrics.

## Permutation strategy

The `perm_strategy` field in `perm_ctrl` determines how each permutation
pass is executed. The two strategies share the same downstream pipeline
(null construction, p-value scoring, FDR correction) — they differ only
in *how* null metric values are produced. Eligible exact engines reuse
preparation across permutations for built-in blocked CV or explicit
custom splits. Other CV specifications retain their per-permutation fold
generation.

- **`"iterate"`** (default):

  Evaluates a **subsampled** set of centers using a prepared engine when
  eligible, or
  [`mvpa_iterate`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
  otherwise. The `subsample` parameter in `perm_ctrl` controls how many
  centers are evaluated per permutation.

  **Null pool size**: `n_perm * n_subsampled_centers`.

  **Best for**: slow classifiers, large brains, limited compute.
  Subsampling reduces the number of centers scored per permutation.

- **`"searchlight"`**:

  Evaluates the **full brain** using a prepared engine when eligible, or
  [`run_searchlight`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  otherwise. The latter retains standard engine dispatch and
  user-defined `run_searchlight.<class>` methods.

  Since the full brain is computed anyway, **all** centers contribute to
  the null distribution (the `subsample` parameter is ignored and a note
  is logged).

  **Null pool size**: `n_perm * all_centers`.

  **Best for**: models with a fast searchlight engine, or when you want
  the richest possible null distribution.

## Examples

``` r
# \donttest{
  ds    <- gen_sample_dataset(c(5, 5, 5), 20, blocks = 2, nlevels = 2)
  cval  <- blocked_cross_validation(ds$design$block_var)
  mdl   <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification",
                      crossval = cval)

  # Strategy 1: subsampled iterator (default, universal)
  pc1   <- permutation_control(n_perm = 10, subsample = 0.2, seed = 1L)
  res1  <- run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc1)
#> INFO [2026-10-07 05:02:28] Running observed searchlight (radius = 3 mm) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-10-07 05:02:28] searchlight engine: sda_fast
#> INFO [2026-10-07 05:02:28] Building searchlight iterator ...
#> INFO [2026-10-07 05:02:28] Subsampling searchlight centers ...
#> INFO [2026-10-07 05:02:28] Using 25 / 125 centers for permutation runs.
#> INFO [2026-10-07 05:02:28] Permutations use the exact 'sda_fast' engine (prepared once).
#> INFO [2026-10-07 05:02:28] Permutation 1 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 2 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 3 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 4 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 5 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 6 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 7 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 8 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 9 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Permutation 10 / 10 (strategy: iterate) ...
#> INFO [2026-10-07 05:02:28] Running null diagnostics for 'Accuracy' ...
#> Null Distribution Diagnostics (n_perm = 10 )
#> -------------------------------------------------- 
#>   nfeatures            [FLAGGED]
#>     rho=-0.187  p=0.0030
#>     Null correlates with nfeatures (p < 0.01); covariate adjustment recommended.
#> -------------------------------------------------- 
#> INFO [2026-10-07 05:02:28] Building adjusted null distribution for 'Accuracy' (adjusted, 5 bins) ...
#> INFO [2026-10-07 05:02:28] Computing p-values for 125 centers ...
#> INFO [2026-10-07 05:02:28] Building p-value spatial maps ...
#> INFO [2026-10-07 05:02:28] Done 'Accuracy'. 0 centers significant at FDR < 0.05 (fdr).

  # Strategy 2: full-brain via run_searchlight (engine-aware)
  pc2   <- permutation_control(n_perm = 5, perm_strategy = "searchlight",
                               seed = 1L)
  res2  <- run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc2)
#> INFO [2026-10-07 05:02:28] Running observed searchlight (radius = 3 mm) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-10-07 05:02:29] searchlight engine: sda_fast
#> INFO [2026-10-07 05:02:29] Building searchlight iterator ...
#> INFO [2026-10-07 05:02:29] Strategy = 'searchlight': full brain computed per permutation. 'subsample' parameter ignored; all 125 centers contribute to null.
#> INFO [2026-10-07 05:02:29] Permutations use the exact 'sda_fast' engine (prepared once).
#> INFO [2026-10-07 05:02:29] Permutation 1 / 5 (strategy: searchlight) ...
#> INFO [2026-10-07 05:02:29] Permutation 2 / 5 (strategy: searchlight) ...
#> INFO [2026-10-07 05:02:29] Permutation 3 / 5 (strategy: searchlight) ...
#> INFO [2026-10-07 05:02:29] Permutation 4 / 5 (strategy: searchlight) ...
#> INFO [2026-10-07 05:02:29] Permutation 5 / 5 (strategy: searchlight) ...
#> INFO [2026-10-07 05:02:29] Running null diagnostics for 'Accuracy' ...
#> Null Distribution Diagnostics (n_perm = 5 )
#> -------------------------------------------------- 
#>   nfeatures            [FLAGGED]
#>     rho=-0.169  p=0.0000
#>     Null correlates with nfeatures (p < 0.01); covariate adjustment recommended.
#> -------------------------------------------------- 
#> INFO [2026-10-07 05:02:29] Building adjusted null distribution for 'Accuracy' (adjusted, 5 bins) ...
#> INFO [2026-10-07 05:02:29] Computing p-values for 125 centers ...
#> INFO [2026-10-07 05:02:29] Building p-value spatial maps ...
#> INFO [2026-10-07 05:02:29] Done 'Accuracy'. 0 centers significant at FDR < 0.05 (fdr).
# }
```
