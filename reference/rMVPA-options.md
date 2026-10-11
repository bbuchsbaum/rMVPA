# rMVPA runtime options

Options read with [`getOption()`](https://rdrr.io/r/base/options.html)
by the package at run time. Set them with `options(name = value)`. Each
default below is used when the option is unset. Byte budgets are plain
numbers of bytes.

## Details

- `rMVPA.regional_mem_budget`:

  Default `2e9` (2 GB). Peak extraction-memory budget, in bytes, for one
  batch of regional ROIs. The batch size is reduced when the estimated
  per-batch memory exceeds it.

- `rMVPA.searchlight_mem_budget`:

  Default `512 * 1024^2` (512 MiB). Memory budget, in bytes, used to
  pick the automatic batch size for searchlight centers.

- `rMVPA.searchlight_min_chunk`:

  Default `16L`. Lower bound, in centers, on the chunk size used when
  searchlight centers are split across workers. A chunk is still capped
  by the centers per worker.

- `rMVPA.gc_each_batch`:

  Default `FALSE`. When `TRUE`, the iterator calls
  [`gc()`](https://rdrr.io/r/base/gc.html) after every batch. Off by
  default because a full collection is costly; use it to lower peak
  memory in long runs.

- `rMVPA.shard_gc_each_batch`:

  Default `FALSE`. Older alias of `rMVPA.gc_each_batch`, honoured by the
  iterator and the shard backend.

- `rMVPA.profile_searchlight`:

  Default `FALSE`. When `TRUE`, searchlight runs record setup and timing
  information in the `"timing"` attribute of the result.

- `rMVPA.debug`:

  Default `FALSE`. When `TRUE`, the searchlight combiner resolution
  reports extra diagnostic information.

- `rMVPA.rsa_perm_cache_bytes`:

  Default `512 * 1024^2` (512 MiB). Maximum size, in bytes, of the RSA
  permutation cache. When the cache for a run would exceed this, it is
  not used.

- `rMVPA.rsa_fast_parallel`:

  Default `TRUE`. Allows the fast RSA searchlight engine to run its
  chunks in parallel when more than one worker is available and there
  are at least 256 centers. Set to `FALSE` to force serial execution.

- `rMVPA.feature_rsa_block_target_bytes`:

  Default `128 * 1024^2` (128 MiB). Target size, in bytes, of one
  feature-RSA block. Non-finite or non-positive values fall back to the
  default.

- `rMVPA.feature_rsa_block_max_cols`:

  Default `64L`. Upper bound on the number of columns in one feature-RSA
  block.

- `rMVPA.feature_rsa_metric_block_rows`:

  Default `128L`. Number of rows processed per block when computing
  feature-RSA metrics. It is clamped to the number of observations.

- `rMVPA.feature_rsa_component_target_bytes`:

  Default `32 * 1024^2` (32 MiB). Target size, in bytes, of one block of
  components when predicting feature-RSA models with a PLS/PCR fit.

Options that are not listed here are not read by rMVPA.
