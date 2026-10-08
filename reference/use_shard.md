# Enable Shared-Memory Backend for an MVPA Model Spec

Prepares the dataset for shared-memory access and tags the model
specification so that
[`mvpa_iterate`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
uses the shard backend instead of the default furrr pipeline.

## Usage

``` r
use_shard(mod_spec)
```

## Arguments

- mod_spec:

  A model specification created by
  [`mvpa_model`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md)
  (or any constructor that produces an object inheriting from
  `"model_spec"`).

## Value

A copy of `mod_spec` with class `"shard_model_spec"` prepended and a
`$shard_data` list attached.

## Details

This is an **experimental** feature. It requires the shard package and a
platform that supports POSIX shared memory (`shm_open`).

## Examples

``` r
if (requireNamespace("shard", quietly = TRUE)) {
  ds <- gen_sample_dataset(c(4, 4, 4), 24, nlevels = 2, blocks = 3)
  cval <- blocked_cross_validation(ds$design$block_var)
  mspec <- mvpa_model(load_model("corclass"), ds$dataset, ds$design,
                      "classification", crossval = cval)
  mspec <- use_shard(mspec)
  results <- run_searchlight(mspec, radius = 2, method = "randomized",
                             niter = 1)
  shard_cleanup(mspec$shard_data)
}
#> INFO [2026-10-08 02:06:45] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-10-08 02:06:45] shard backend [volumetric]: shared 24 x 64 matrix (64 masked voxels)
#> INFO [2026-10-08 02:06:45] searchlight engine: general-purpose iterator (no eligible fast path)
#> INFO [2026-10-08 02:06:45] Running randomized searchlight with radius = 2 and niter = 1
#> INFO [2026-10-08 02:06:45] Starting randomized searchlight analysis:
#> INFO [2026-10-08 02:06:45] - Radius: 2
#> INFO [2026-10-08 02:06:45] - Iterations: 1
#> INFO [2026-10-08 02:06:45] 
#> Iteration 1/1
#> INFO [2026-10-08 02:06:45] Using automatic searchlight batch size 8 for 8 centers (memory budget 512.0 MiB).
#> INFO [2026-10-08 02:06:45] 
#> MVPA Iteration Complete
#> - Total ROIs: 8
#> - Processed: 8
#> - Skipped: 0
#> INFO [2026-10-08 02:06:45] searchlight (randomized): 8 ROIs processed (success=7, errors=1)
#> WARN [2026-10-08 02:06:45] searchlight (randomized): 1 of 8 ROIs failed (12.5%)
#> WARN [2026-10-08 02:06:45]   - [1 ROIs] error: less than 2 features
```
