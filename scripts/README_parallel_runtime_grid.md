# Parallel Runtime Sweep

`scripts/sweep_parallel_runtime_grid.R` is a small HPC triage harness for
finding bad interactions between:

- `future` backends (`sequential`, `multisession`, `multicore`, `callr`,
  `mirai_multisession`)
- worker counts
- OpenMP thread counts
- BLAS thread counts
- `rMVPA` batch sizes

The parent process launches each configuration in a fresh `Rscript` child,
applies a wall-clock timeout, captures per-run logs, and writes both raw and
summary CSV outputs. This is intended for debugging hangs where jobs appear
busy but never complete.

## Why a fresh child per configuration?

Thread-related env vars such as `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and
`MKL_NUM_THREADS` need to be set before native libraries initialize. Running
each combination in a fresh `Rscript` process makes those settings meaningful
and keeps one bad configuration from poisoning the rest of the sweep.

## Quick Start

From the package root:

```bash
env \
  RMVPA_HPC_SWEEP_DRY_RUN=true \
  Rscript scripts/sweep_parallel_runtime_grid.R
```

This writes a manifest only. It does not run any analyses.

For a small real sweep:

```bash
env \
  RMVPA_HPC_SWEEP_ANALYSES=regional \
  RMVPA_HPC_SWEEP_BACKENDS=sequential,multisession,multicore \
  RMVPA_HPC_SWEEP_WORKER_COUNTS=1,2,4 \
  RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS=1,2,4 \
  RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS=1,2 \
  RMVPA_HPC_SWEEP_BATCH_SIZES=auto,4 \
  RMVPA_HPC_SWEEP_REP=2 \
  RMVPA_HPC_SWEEP_TIMEOUT_SECONDS=300 \
  RMVPA_HPC_SWEEP_OUT=.tmp/hpc_parallel_sweep_summary.csv \
  RMVPA_HPC_SWEEP_OUT_RAW=.tmp/hpc_parallel_sweep_raw.csv \
  RMVPA_HPC_SWEEP_LOG_DIR=.tmp/hpc_parallel_sweep_logs \
  Rscript scripts/sweep_parallel_runtime_grid.R
```

## Main Environment Variables

- `RMVPA_HPC_SWEEP_ANALYSES`
  Comma-separated analysis types. Supported: `regional`, `searchlight`.
- `RMVPA_HPC_SWEEP_BACKENDS`
  Comma-separated future backends. Supported:
  `sequential,multisession,multicore,callr,mirai_multisession`.
- `RMVPA_HPC_SWEEP_WORKER_COUNTS`
  Comma-separated worker counts, e.g. `1,2,4,8`.
- `RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS`
  Comma-separated `OMP_NUM_THREADS` values.
- `RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS`
  Comma-separated values applied to:
  `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `BLIS_NUM_THREADS`,
  `VECLIB_MAXIMUM_THREADS`.
- `RMVPA_HPC_SWEEP_BATCH_SIZES`
  Comma-separated batch sizes. Use `auto` and/or positive integers.
- `RMVPA_HPC_SWEEP_REP`
  Number of repetitions per configuration.
- `RMVPA_HPC_SWEEP_TIMEOUT_SECONDS`
  Hard wall-clock timeout per child run.
- `RMVPA_HPC_SWEEP_RADIUS`
  Searchlight radius for `searchlight` sweeps. Ignored for `regional`.
- `RMVPA_HPC_SWEEP_OUT`
  Summary CSV path.
- `RMVPA_HPC_SWEEP_OUT_RAW`
  Raw per-run CSV path.
- `RMVPA_HPC_SWEEP_LOG_DIR`
  Directory for per-run logs and child result files.
- `RMVPA_HPC_SWEEP_DRY_RUN`
  If `true`, only write the planned manifest.

## Output Files

- Summary CSV (`RMVPA_HPC_SWEEP_OUT`)
  One row per configuration, aggregating repeated runs.
- Raw CSV (`RMVPA_HPC_SWEEP_OUT_RAW`)
  One row per actual run.
- Log directory (`RMVPA_HPC_SWEEP_LOG_DIR`)
  One `.log` file per run, useful when a configuration errors or times out.

Raw rows include:

- `status`: `ok`, `error`, `timeout`, `crash`, or `skip`
- `elapsed_seconds`
- `future_backend`
- `workers`
- `omp_threads`
- `blas_threads`
- `batch_size`
- `message`
- `log_file`

## Interpreting Failures

- `skip`
  The backend was unavailable in the current R installation. Example:
  `future.mirai` not installed.
- `error`
  The child process started and returned an R error.
- `timeout`
  The child stayed alive longer than `RMVPA_HPC_SWEEP_TIMEOUT_SECONDS` and was
  terminated.
- `crash`
  The child exited without producing a normal result row. This is a useful
  signal for low-level runtime problems.

For `error`, `timeout`, and `crash`, inspect the corresponding `log_file`.

## Suggested HPC Workflow

1. Start with a dry run to confirm the grid size and output locations.
2. Run a very small sweep first:
   `workers=1,2`, `OMP=1,2`, `BLAS=1`, `REP=1`.
3. Add backends gradually:
   `sequential` -> `multisession` -> `multicore` -> `callr` -> `mirai`.
4. Increase thread counts only after the low-thread grid is stable.
5. Treat `workers * OMP threads * BLAS threads` as a potential oversubscription
   budget, not independent knobs.

## Example SLURM Job

```bash
#!/bin/bash
#SBATCH --job-name=rmvpa-sweep
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/rmvpa-sweep-%j.out

module load R
cd /path/to/rMVPA

mkdir -p .tmp/hpc_logs logs

env \
  RMVPA_HPC_SWEEP_ANALYSES=regional \
  RMVPA_HPC_SWEEP_BACKENDS=sequential,multisession,multicore \
  RMVPA_HPC_SWEEP_WORKER_COUNTS=1,2,4,8 \
  RMVPA_HPC_SWEEP_OMP_THREAD_COUNTS=1,2,4 \
  RMVPA_HPC_SWEEP_BLAS_THREAD_COUNTS=1,2,4 \
  RMVPA_HPC_SWEEP_BATCH_SIZES=auto,4 \
  RMVPA_HPC_SWEEP_REP=1 \
  RMVPA_HPC_SWEEP_TIMEOUT_SECONDS=300 \
  RMVPA_HPC_SWEEP_OUT=.tmp/hpc_parallel_sweep_summary.csv \
  RMVPA_HPC_SWEEP_OUT_RAW=.tmp/hpc_parallel_sweep_raw.csv \
  RMVPA_HPC_SWEEP_LOG_DIR=.tmp/hpc_logs \
  Rscript scripts/sweep_parallel_runtime_grid.R
```

## Notes

- `mirai_multisession` is marked as `skip` unless `future.mirai` is installed.
- `callr` is marked as `skip` unless `future.callr` is installed.
- `multicore` is marked as `skip` when `future::supportsMulticore()` is false.
- The script currently targets Unix-like systems. It is not designed for
  Windows.
