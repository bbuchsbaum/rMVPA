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
in *how* null metric values are produced. Neither strategy contains any
engine-specific branching.

- **`"iterate"`** (default):

  Calls
  [`mvpa_iterate`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
  on a **subsampled** set of centers. This is the generic per-ROI
  iterator that works with every model type and every dataset class. The
  `subsample` parameter in `perm_ctrl` controls how many centers are
  evaluated per permutation.

  **Null pool size**: `n_perm * n_subsampled_centers`.

  **Best for**: slow classifiers, large brains, limited compute. The
  subsampling gives a 5–20\\\times\\ speedup over a full-brain pass.

- **`"searchlight"`**:

  Calls
  [`run_searchlight`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  on the **full brain** for each permutation. Because the call goes
  through the standard `run_searchlight` S3 dispatch, it automatically
  benefits from any fast engine the model qualifies for (e.g.\\ SWIFT,
  dual-LDA) *and* from any user-defined `run_searchlight.<class>`
  method. No engine-specific code exists here — it is purely the
  standard `run_searchlight` call.

  Since the full brain is computed anyway, **all** centers contribute to
  the null distribution (the `subsample` parameter is ignored and a note
  is logged). This yields a richer null and therefore better-calibrated
  p-values.

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
#> INFO [2026-09-05 00:59:41] Running observed searchlight (radius = 3 mm) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-09-05 00:59:41] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-09-05 00:59:41] Running standard searchlight with radius = 3
#> INFO [2026-09-05 00:59:41] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-09-05 00:59:41] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-09-05 00:59:41] creating standard searchlight
#> INFO [2026-09-05 00:59:41] running standard searchlight iterator
#> INFO [2026-09-05 00:59:41] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:44] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-09-05 00:59:45] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-09-05 00:59:45] Building searchlight iterator ...
#> INFO [2026-09-05 00:59:45] Subsampling searchlight centers ...
#> INFO [2026-09-05 00:59:45] Using 25 / 125 centers for permutation runs.
#> INFO [2026-09-05 00:59:45] Permutation 1 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:45] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:47] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 00:59:47] Permutation 2 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:47] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:49] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 00:59:49] Permutation 3 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:49] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:51] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 00:59:51] Permutation 4 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:51] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:53] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 00:59:53] Permutation 5 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:53] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:55] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 00:59:55] Permutation 6 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:55] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:57] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 00:59:57] Permutation 7 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:57] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 00:59:59] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 00:59:59] Permutation 8 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 00:59:59] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:01] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 01:00:01] Permutation 9 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 01:00:01] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:03] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 01:00:03] Permutation 10 / 10 (strategy: iterate) ...
#> INFO [2026-09-05 01:00:03] Using automatic searchlight batch size 25 for 25 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:05] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-09-05 01:00:06] Running null diagnostics for 'Accuracy' ...
#> Null Distribution Diagnostics (n_perm = 10 )
#> -------------------------------------------------- 
#>   nfeatures            [OK]
#>     rho=-0.063  p=0.3224
#>     No significant correlation with nfeatures.
#> -------------------------------------------------- 
#> INFO [2026-09-05 01:00:06] Building adjusted null distribution for 'Accuracy' (adjusted, 5 bins) ...
#> INFO [2026-09-05 01:00:06] Computing p-values for 125 centers ...
#> INFO [2026-09-05 01:00:06] Building p-value spatial maps ...
#> INFO [2026-09-05 01:00:06] Done 'Accuracy'. 0 centers significant at FDR < 0.05 (fdr).

  # Strategy 2: full-brain via run_searchlight (engine-aware)
  pc2   <- permutation_control(n_perm = 5, perm_strategy = "searchlight",
                               seed = 1L)
  res2  <- run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc2)
#> INFO [2026-09-05 01:00:06] Running observed searchlight (radius = 3 mm) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-09-05 01:00:06] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-09-05 01:00:06] Running standard searchlight with radius = 3
#> INFO [2026-09-05 01:00:06] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-09-05 01:00:06] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-09-05 01:00:06] creating standard searchlight
#> INFO [2026-09-05 01:00:06] running standard searchlight iterator
#> INFO [2026-09-05 01:00:06] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:09] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-09-05 01:00:10] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-09-05 01:00:10] Building searchlight iterator ...
#> INFO [2026-09-05 01:00:10] Strategy = 'searchlight': full brain computed per permutation. 'subsample' parameter ignored; all 125 centers contribute to null.
#> INFO [2026-09-05 01:00:10] Permutation 1 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-09-05 01:00:10] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-09-05 01:00:10] Running standard searchlight with radius = 3
#> INFO [2026-09-05 01:00:10] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-09-05 01:00:10] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-09-05 01:00:10] creating standard searchlight
#> INFO [2026-09-05 01:00:10] running standard searchlight iterator
#> INFO [2026-09-05 01:00:10] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:14] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-09-05 01:00:14] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-09-05 01:00:14] Permutation 2 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-09-05 01:00:14] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-09-05 01:00:14] Running standard searchlight with radius = 3
#> INFO [2026-09-05 01:00:14] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-09-05 01:00:14] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-09-05 01:00:14] creating standard searchlight
#> INFO [2026-09-05 01:00:14] running standard searchlight iterator
#> INFO [2026-09-05 01:00:14] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:18] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-09-05 01:00:18] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-09-05 01:00:18] Permutation 3 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-09-05 01:00:18] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-09-05 01:00:18] Running standard searchlight with radius = 3
#> INFO [2026-09-05 01:00:18] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-09-05 01:00:18] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-09-05 01:00:18] creating standard searchlight
#> INFO [2026-09-05 01:00:18] running standard searchlight iterator
#> INFO [2026-09-05 01:00:18] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:23] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-09-05 01:00:23] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-09-05 01:00:23] Permutation 4 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-09-05 01:00:23] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-09-05 01:00:23] Running standard searchlight with radius = 3
#> INFO [2026-09-05 01:00:23] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-09-05 01:00:23] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-09-05 01:00:23] creating standard searchlight
#> INFO [2026-09-05 01:00:23] running standard searchlight iterator
#> INFO [2026-09-05 01:00:23] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:27] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-09-05 01:00:27] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-09-05 01:00:27] Permutation 5 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-09-05 01:00:27] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-09-05 01:00:27] Running standard searchlight with radius = 3
#> INFO [2026-09-05 01:00:27] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-09-05 01:00:27] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-09-05 01:00:27] creating standard searchlight
#> INFO [2026-09-05 01:00:27] running standard searchlight iterator
#> INFO [2026-09-05 01:00:27] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-09-05 01:00:31] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-09-05 01:00:32] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-09-05 01:00:32] Running null diagnostics for 'Accuracy' ...
#> Null Distribution Diagnostics (n_perm = 5 )
#> -------------------------------------------------- 
#>   nfeatures            [OK]
#>     rho=-0.084  p=0.0350
#>     No significant correlation with nfeatures.
#> -------------------------------------------------- 
#> INFO [2026-09-05 01:00:32] Building adjusted null distribution for 'Accuracy' (adjusted, 5 bins) ...
#> INFO [2026-09-05 01:00:32] Computing p-values for 125 centers ...
#> INFO [2026-09-05 01:00:32] Building p-value spatial maps ...
#> INFO [2026-09-05 01:00:32] Done 'Accuracy'. 0 centers significant at FDR < 0.05 (fdr).
# }
```
