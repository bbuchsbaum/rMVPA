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

  An `mvpa_model` or compatible model specification.

- observed:

  Optional pre-computed searchlight result (output of
  [`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)).
  If `NULL`, it is computed internally. Can also be a named numeric
  vector of observed metric values indexed by center ID.

- radius:

  Numeric. Searchlight radius in mm (default 8).

- method:

  Character. Searchlight method passed to
  [`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  when `observed` is `NULL` or when `perm_strategy = "searchlight"`.

- perm_ctrl:

  A
  [`permutation_control`](http://bbuchsbaum.github.io/rMVPA/reference/permutation_control.md)
  object.

- metric:

  Character. Which performance metric to use for inference. If `NULL`
  (default), the first metric is used.

- ...:

  Additional arguments forwarded to
  [`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  (observed pass and `"searchlight"` permutations). For `"iterate"`
  permutations, only arguments that are formal parameters of
  [`mvpa_iterate`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
  are forwarded; other keys are ignored for that path to avoid
  argument-mismatch failures (e.g., `engine = "legacy"` is meaningful
  for `run_searchlight` but not for `mvpa_iterate`).

## Value

A `permutation_result` S3 object containing:

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

## Permutation strategy

The `perm_strategy` field in `perm_ctrl` determines how each permutation
pass is executed. The two strategies share the same downstream pipeline
(null construction, p-value scoring, FDR correction) — they differ only
in *how* null metric values are produced. Neither strategy contains any
engine-specific branching.

- **`"iterate"`** (default):

  Calls
  [`mvpa_iterate`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
  on a **subsampled** set of centers. This is the generic per-ROI
  iterator that works with every model type and every dataset class. The
  `subsample` parameter in `perm_ctrl` controls how many centers are
  evaluated per permutation.

  **Null pool size**: `n_perm * n_subsampled_centers`.

  **Best for**: slow classifiers, large brains, limited compute. The
  subsampling gives a 5–20\\\times\\ speedup over a full-brain pass.

- **`"searchlight"`**:

  Calls
  [`run_searchlight`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
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
#> INFO [2026-04-11 19:48:07] Running observed searchlight (radius = 3 mm) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-04-11 19:48:07] searchlight engine: legacy (no eligible fast path)
#> INFO [2026-04-11 19:48:07] Running standard searchlight with radius = 3
#> INFO [2026-04-11 19:48:07] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-04-11 19:48:07] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-04-11 19:48:07] creating standard searchlight
#> INFO [2026-04-11 19:48:07] running standard searchlight iterator
#> INFO [2026-04-11 19:48:18] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-04-11 19:48:18] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-04-11 19:48:18] Building searchlight iterator ...
#> INFO [2026-04-11 19:48:18] Subsampling searchlight centers ...
#> INFO [2026-04-11 19:48:18] Using 25 / 125 centers for permutation runs.
#> INFO [2026-04-11 19:48:18] Permutation 1 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:20] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:20] Permutation 2 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:23] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:23] Permutation 3 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:25] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:25] Permutation 4 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:28] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:28] Permutation 5 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:30] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:30] Permutation 6 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:33] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:33] Permutation 7 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:35] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:35] Permutation 8 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:38] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:38] Permutation 9 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:40] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:40] Permutation 10 / 10 (strategy: iterate) ...
#> INFO [2026-04-11 19:48:43] 
#> MVPA Iteration Complete
#> - Total ROIs: 25
#> - Processed: 25
#> - Skipped: 0
#> INFO [2026-04-11 19:48:43] Running null diagnostics ...
#> Null Distribution Diagnostics (n_perm = 10 )
#> -------------------------------------------------- 
#>   nfeatures            [OK]
#>     rho=-0.023  p=0.7178
#>     No significant correlation with nfeatures.
#> -------------------------------------------------- 
#> INFO [2026-04-11 19:48:43] Building adjusted null distribution (adjusted, 5 bins) ...
#> INFO [2026-04-11 19:48:43] Computing p-values for 125 centers ...
#> INFO [2026-04-11 19:48:43] Building p-value spatial maps ...
#> INFO [2026-04-11 19:48:43] Done. 0 centers significant at FDR < 0.05 (fdr).

  # Strategy 2: full-brain via run_searchlight (engine-aware)
  pc2   <- permutation_control(n_perm = 5, perm_strategy = "searchlight",
                               seed = 1L)
  res2  <- run_permutation_searchlight(mspec, radius = 3, perm_ctrl = pc2)
#> INFO [2026-04-11 19:48:43] Running observed searchlight (radius = 3 mm) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-04-11 19:48:43] searchlight engine: legacy (no eligible fast path)
#> INFO [2026-04-11 19:48:43] Running standard searchlight with radius = 3
#> INFO [2026-04-11 19:48:43] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-04-11 19:48:43] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-04-11 19:48:43] creating standard searchlight
#> INFO [2026-04-11 19:48:43] running standard searchlight iterator
#> INFO [2026-04-11 19:48:53] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-04-11 19:48:54] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-04-11 19:48:54] Building searchlight iterator ...
#> INFO [2026-04-11 19:48:54] Strategy = 'searchlight': full brain computed per permutation. 'subsample' parameter ignored; all 125 centers contribute to null.
#> INFO [2026-04-11 19:48:54] Permutation 1 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-04-11 19:48:54] searchlight engine: legacy (no eligible fast path)
#> INFO [2026-04-11 19:48:54] Running standard searchlight with radius = 3
#> INFO [2026-04-11 19:48:54] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-04-11 19:48:54] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-04-11 19:48:54] creating standard searchlight
#> INFO [2026-04-11 19:48:54] running standard searchlight iterator
#> INFO [2026-04-11 19:49:04] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-04-11 19:49:05] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-04-11 19:49:05] Permutation 2 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-04-11 19:49:05] searchlight engine: legacy (no eligible fast path)
#> INFO [2026-04-11 19:49:05] Running standard searchlight with radius = 3
#> INFO [2026-04-11 19:49:05] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-04-11 19:49:05] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-04-11 19:49:05] creating standard searchlight
#> INFO [2026-04-11 19:49:05] running standard searchlight iterator
#> INFO [2026-04-11 19:49:15] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-04-11 19:49:16] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-04-11 19:49:16] Permutation 3 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-04-11 19:49:16] searchlight engine: legacy (no eligible fast path)
#> INFO [2026-04-11 19:49:16] Running standard searchlight with radius = 3
#> INFO [2026-04-11 19:49:16] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-04-11 19:49:16] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-04-11 19:49:16] creating standard searchlight
#> INFO [2026-04-11 19:49:16] running standard searchlight iterator
#> INFO [2026-04-11 19:49:26] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-04-11 19:49:27] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-04-11 19:49:27] Permutation 4 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-04-11 19:49:27] searchlight engine: legacy (no eligible fast path)
#> INFO [2026-04-11 19:49:27] Running standard searchlight with radius = 3
#> INFO [2026-04-11 19:49:27] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-04-11 19:49:27] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-04-11 19:49:27] creating standard searchlight
#> INFO [2026-04-11 19:49:27] running standard searchlight iterator
#> INFO [2026-04-11 19:49:37] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-04-11 19:49:38] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-04-11 19:49:38] Permutation 5 / 5 (strategy: searchlight) ...
#> Warning: run_searchlight preflight reported 0 failure(s) and 1 warning(s).
#> block_count: Only 2 blocks found. Leave-one-block-out CV will have only 2 folds, providing limited evaluation with high variance estimates.
#> INFO [2026-04-11 19:49:38] searchlight engine: legacy (no eligible fast path)
#> INFO [2026-04-11 19:49:38] Running standard searchlight with radius = 3
#> INFO [2026-04-11 19:49:38] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-04-11 19:49:38] shard backend [volumetric]: shared 20 x 125 matrix (125 masked voxels)
#> INFO [2026-04-11 19:49:38] creating standard searchlight
#> INFO [2026-04-11 19:49:38] running standard searchlight iterator
#> INFO [2026-04-11 19:49:48] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-04-11 19:49:48] searchlight (standard): 125 ROIs processed (success=125, errors=0)
#> INFO [2026-04-11 19:49:48] Running null diagnostics ...
#> Null Distribution Diagnostics (n_perm = 5 )
#> -------------------------------------------------- 
#>   nfeatures            [OK]
#>     rho=-0.077  p=0.0545
#>     No significant correlation with nfeatures.
#> -------------------------------------------------- 
#> INFO [2026-04-11 19:49:48] Building adjusted null distribution (adjusted, 5 bins) ...
#> INFO [2026-04-11 19:49:48] Computing p-values for 125 centers ...
#> INFO [2026-04-11 19:49:48] Building p-value spatial maps ...
#> INFO [2026-04-11 19:49:48] Done. 0 centers significant at FDR < 0.05 (fdr).
# }
```
