# Experimental Shared-Memory Backend for Parallel MVPA

Uses the shard package to place the full data matrix in POSIX shared
memory, enabling zero-copy parallel ROI extraction across furrr workers.
Activated per model spec via
[`use_shard`](http://bbuchsbaum.github.io/rMVPA/reference/use_shard.md).

## Details

When enabled, the default serial ROI-extraction loop in
[`mvpa_iterate`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
is skipped. Instead, each worker receives lightweight ALTREP handles
that point into the shared-memory segment and extracts its own ROI
columns on demand. This eliminates the two main memory bottlenecks of
the default pipeline:

1.  Serial `extract_roi()` in the main process.

2.  Serialisation of large ROI matrices to workers via furrr.

Searchlight batches carry only neighborhood indices, so automatic shard
batches are larger than default-backend batches to amortize scheduling
overhead. Full garbage collection is not forced between shard batches by
default; set `options(rMVPA.shard_gc_each_batch = TRUE)` to restore that
conservative behavior when diagnosing memory pressure.

## Supported dataset types

- `mvpa_image_dataset` (volumetric)

- `mvpa_surface_dataset` (cortical surface)

- `mvpa_clustered_dataset` (parcellated)

## Cleanup

Shared-memory segments are released automatically at the end of
[`mvpa_iterate()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md).
You can also call `shard_cleanup(mod_spec$shard_data)` explicitly.
