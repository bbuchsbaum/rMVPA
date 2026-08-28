# Parallel Runtime Sweep

`scripts/sweep_parallel_runtime_grid.R` characterizes three independent layers
of parallel execution:

1. the Future process topology (`sequential`, `multisession`, `multicore`,
   `callr`, or `mirai_multisession`);
2. the rMVPA data-transport backend (`default` or `shard`); and
3. native thread limits for OpenMP and BLAS.

Future and furrr schedule work. The rMVPA backend controls how ROI data reach
that work. Native libraries may create threads inside each Future worker. Do
not treat these as interchangeable settings.

The parent launches every configuration in a fresh `Rscript`, applies a hard
timeout, captures a log, and writes raw and summary CSVs. A successful child
also records a canonical numerical result, a result signature, task-frame
sizes, phase timings, and sampled process-tree RSS.

## Why a fresh child per configuration?

Thread variables such as `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and
`MKL_NUM_THREADS` should be set before a native runtime initializes. A fresh
process makes those settings meaningful and prevents a crash or hang in one
configuration from contaminating later runs.

Fresh processes also expose startup costs. For a short analysis, a socket or
mirai plan can spend more time starting and loading workers than fitting the
model. That is part of the observed run, but it does not predict performance
for a long-lived worker pool.

## Quick start

From the package root, write a manifest without running analyses:

```bash
env \
  RMVPA_HPC_SWEEP_DRY_RUN=true \
  RMVPA_HPC_SWEEP_DATA_BACKENDS=default,shard \
  Rscript scripts/sweep_parallel_runtime_grid.R
```

For a small real sweep:

```bash
env \
  RMVPA_HPC_SWEEP_ANALYSES=regional \
  RMVPA_HPC_SWEEP_BACKENDS=sequential,multisession,multicore \
  RMVPA_HPC_SWEEP_DATA_BACKENDS=default,shard \
  RMVPA_HPC_SWEEP_WORKER_COUNTS=1,2 \
  RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS=1 \
  RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS=1 \
  RMVPA_HPC_SWEEP_BATCH_SIZES=4 \
  RMVPA_HPC_SWEEP_REP=2 \
  RMVPA_HPC_SWEEP_TIMEOUT_SECONDS=300 \
  RMVPA_HPC_SWEEP_OUT=.tmp/parallel_sweep_summary.csv \
  RMVPA_HPC_SWEEP_OUT_RAW=.tmp/parallel_sweep_raw.csv \
  RMVPA_HPC_SWEEP_LOG_DIR=.tmp/parallel_sweep_logs \
  Rscript scripts/sweep_parallel_runtime_grid.R
```

Use `batch_size` small enough to produce multiple batches when the purpose is
to exercise multiple workers. The publication driver refuses a receipt if
`n_batches < 2`.

## Benchmark the exact checkout

`pkgload::load_all()` is useful for debugging a sequential child, but a clean
installed package is the defensible target for process-based measurements.
Build the checkout, install the tarball into a temporary library, and put that
library first:

```bash
R CMD build --no-build-vignettes --no-manual .
benchmark_lib=$(mktemp -d /tmp/rmvpa-benchmark-lib.XXXXXX)
R CMD INSTALL --library="$benchmark_lib" rMVPA_*.tar.gz
R_LIBS="$benchmark_lib" \
  Rscript inst/benchmarks/parallel_runtime_benchmark.R
```

The benchmark driver fixes the seed and fixture, requests five batches,
configures two workers for non-sequential plans, pins OpenMP and BLAS to one
thread, runs three fresh-process repetitions, and refuses to publish unless
all scheduler/data-backend results match to `1e-12` or better.

## Main environment variables

- `RMVPA_HPC_SWEEP_ANALYSES`: `regional` and/or `searchlight`.
- `RMVPA_HPC_SWEEP_BACKENDS`: `sequential`, `multisession`, `multicore`,
  `callr`, and/or `mirai_multisession`.
- `RMVPA_HPC_SWEEP_DATA_BACKENDS`: `default` and/or `shard`.
- `RMVPA_HPC_SWEEP_LOAD_MODE`: `installed` for evidence; `source` for local
  debugging.
- `RMVPA_HPC_SWEEP_WORKER_COUNTS`: requested Future worker counts.
- `RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS`: `OMP_NUM_THREADS` values.
- `RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS`: values applied to
  `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `BLIS_NUM_THREADS`, and
  `VECLIB_MAXIMUM_THREADS`.
- `RMVPA_HPC_SWEEP_BATCH_SIZES`: `auto` and/or positive integers.
- `RMVPA_HPC_SWEEP_REP`: fresh-process repetitions per configuration.
- `RMVPA_HPC_SWEEP_TIMEOUT_SECONDS`: hard timeout per child.
- `RMVPA_HPC_SWEEP_DIMS`: three comma-separated spatial dimensions.
- `RMVPA_HPC_SWEEP_NOBS`, `RMVPA_HPC_SWEEP_NLEVELS`,
  `RMVPA_HPC_SWEEP_BLOCKS`, and `RMVPA_HPC_SWEEP_N_REGIONS`: fixture shape.
- `RMVPA_HPC_SWEEP_MODEL`: registered rMVPA model ID.
- `RMVPA_HPC_SWEEP_RADIUS`: searchlight radius; ignored for regional runs.
- `RMVPA_HPC_SWEEP_OUT`, `RMVPA_HPC_SWEEP_OUT_RAW`, and
  `RMVPA_HPC_SWEEP_LOG_DIR`: output locations.
- `RMVPA_HPC_SWEEP_DRY_RUN`: write the planned grid without analyses.

## Reading the receipt

The raw CSV contains one row per run. Important fields are:

- `status`: `ok`, `error`, `timeout`, `crash`, or `skip`;
- `configured_workers` and `n_batches`: evidence that the intended path ran;
- `result_keys`, `result_values`, and `result_signature`: the correctness
  receipt;
- `analysis_seconds`, `roi_extract_seconds`, and `run_future_seconds`: nested
  timing phases;
- `data_bytes`, `frame_bytes_sum`, and `frame_bytes_max`: R object sizes for
  the source dataset and serialized task frames; and
- `peak_tree_rss_bytes`, `peak_process_count`, and `memory_samples`: sampled
  process accounting.

Task-frame bytes are the robust structural comparison. The default backend
strips the full dataset from the worker model specification, but the controller
extracts ROI matrices and places them in task frames. Shard instead shares a
matrix and sends lightweight indices. Thus, “multisession copies the entire
dataset” is not an accurate description of rMVPA's default path, while
“multisession has independent R processes and receives serialized task data”
is accurate.

Do not interpret summed RSS as unique physical memory. It can double-count
shared and copy-on-write pages, and detached worker implementations may not
remain descendants of the controller. The checked-in publication receipt does
not use RSS for its memory claim.

## Interpreting failures

- `skip`: the backend is unavailable. For example, `future.mirai` or `shard`
  is not installed, POSIX shared memory is unavailable, or multicore is not
  supported.
- `error`: the child returned an R error.
- `timeout`: the child stayed alive beyond the configured deadline and was
  terminated.
- `crash`: the child exited without producing a result row.

Inspect the corresponding `log_file` for every non-`ok` row. Matching status
is not enough: compare canonical results within each fixture repetition.

## OpenMP and oversubscription

`OMP_NUM_THREADS` controls native OpenMP regions; it is not a Future worker
count. A conservative upper budget is roughly Future workers times the maximum
native threads used by each worker, with BLAS or other thread pools accounted
for separately. The product is a capacity warning, not a promise that all
pools are active simultaneously.

An rMVPA workload that never enters an OpenMP region cannot validate fork plus
OpenMP safety merely by sweeping `OMP_NUM_THREADS`. Use the sweep to isolate
and reproduce a known native workload, not to certify a forked process after a
few successful runs. R Core strongly discourages forked processing with
multithreaded libraries; prefer separate R sessions and pin native threads
unless the complete native stack is known to be fork-safe.

## Suggested HPC workflow

1. Use `parallelly::availableCores()` and the scheduler allocation rather than
   the machine's physical core count.
2. Start with `multisession`, `default`, two workers, and all native pools at
   one thread.
3. Confirm more than one batch, finite outputs, and numerical parity with the
   sequential result.
4. Add `shard` as a separate data-transport comparison.
5. Add `mirai_multisession` only after its package and worker lifecycle are
   verified on the target cluster.
6. Treat `multicore` as an explicit Unix-only experiment, not the safe default.
7. Scale workers, batch size, and native threads one axis at a time.

`future.mirai::mirai_multisession` still uses separate R processes. It changes
scheduling and communication; it does not itself turn rMVPA data into shared
memory. It can be combined with the shard backend because these are orthogonal
layers.
