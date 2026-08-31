# Choose a Safe Parallel Runtime for rMVPA

A large regional or searchlight analysis can be parallel in three
different ways at once. Future decides **which R processes run tasks**.
rMVPA decides **how image data reach those tasks**. A fitted model or
one of its dependencies may also create **native OpenMP or BLAS threads
inside each R process**.

These controls solve different problems. `future::plan(multisession)`
does not activate shard; `backend = "shard"` does not select a Future
plan; and `OMP_NUM_THREADS` is not an rMVPA worker count. Most
surprising behavior comes from changing two or three layers without
accounting for their product.

For a conservative local or cluster baseline, use separate R sessions,
make the worker count respect the current allocation, and pin native
thread pools to one before starting R:

``` bash
OMP_NUM_THREADS=1 \
OPENBLAS_NUM_THREADS=1 \
MKL_NUM_THREADS=1 \
VECLIB_MAXIMUM_THREADS=1 \
R
```

``` r

run_safely <- function(model_spec, radius = 3) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  workers <- max(1L, parallelly::availableCores(omit = 1L))
  future::plan(future::multisession, workers = workers)

  run_searchlight(
    model_spec,
    radius = radius,
    backend = "default"
  )
}
```

That is the portable starting point, not a claim that it is fastest for
every analysis. Measure the real model, region sizes, and batch
structure on the machine where the analysis will run.

## Keep the three layers separate

| layer | controlled_by | examples | main_cost |
|:---|:---|:---|:---|
| Process topology | future::plan() | sequential, multisession, mirai, multicore | process startup, scheduling, R-session memory |
| Data transport | backend in run_regional()/run_searchlight() | default, shard, auto | ROI extraction, serialization, shared-memory setup |
| Native threads | OMP/BLAS runtime settings | OMP_NUM_THREADS, BLAS-specific limits | oversubscription and fork safety |

Three independent runtime layers in rMVPA. {.table}

furrr is the mapping layer used by rMVPA’s general iterator. It follows
the active Future plan; it is not another process backend. Batch size
matters because a batch is prepared by the controller and then mapped
through furrr. If an analysis creates only one batch, requesting eight
workers does not make eight batches appear.

## Which Future plan should you choose?

| plan | execution | use_when | principal_cost |
|:---|:---|:---|:---|
| sequential | current R process | debugging, a numerical reference, or one allocated core | no process parallelism |
| multisession | persistent independent R sessions | portable local/HPC parallelism and the safest default | startup, serialization, and per-session state |
| mirai_multisession | independent mirai daemon processes | mirai scheduling is already deployed and measured | optional dependency, startup, transport, and daemon lifecycle |
| multicore | forked copies of the current process | Unix-only stack is known to be fork-safe | unsafe with some GUIs and multithreaded native runtimes |
| callr | fresh callr process per future | strong task isolation or prompt worker memory reclamation matters | high startup cost for short futures |

Future process topologies relevant to rMVPA. {.table}

`multisession` is the default recommendation because its workers are
separate R sessions and it is available across platforms. Those sessions
load their own R runtime, packages, model state, and task data, so their
memory is not free. Worker startup can dominate small jobs; reuse a plan
for a sequence of analyses when that is scientifically and operationally
appropriate.

[`future.mirai::mirai_multisession`](https://future.mirai.futureverse.org/reference/mirai_multisession.html)
also uses separate R processes. It changes the scheduler and
communication mechanism; it is not shared memory. It can be combined
with `backend = "shard"` because the two settings operate at different
layers:

``` r

run_with_mirai_and_shard <- function(model_spec, region_mask) {
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)

  future::plan(future.mirai::mirai_multisession, workers = 4)
  run_regional(
    model_spec,
    region_mask,
    backend = "shard"
  )
}
```

`multicore` can start quickly and initially benefits from copy-on-write
pages, but it is not available on Windows and is disabled in some front
ends. More importantly, R Core strongly discourages forked processing
when the parent has used multithreaded libraries: locks and runtime
state can be inherited in an inconsistent form, leading to hangs or
crashes. Copy-on-write is a memory mechanism, not a thread-safety
guarantee. Use multicore only as an explicit, Unix-only experiment after
auditing the entire native stack.

## What do default and shard actually move?

The default backend does **not** export the complete rMVPA dataset with
every future. Before furrr sees the model specification, rMVPA removes
its dataset. The controller extracts each requested ROI, however, and
the resulting ROI matrices live in the task frame sent to workers. Large
or overlapping ROIs can therefore create substantial serialization and
transient task memory even though the original dataset was stripped.

Shard changes that transport path. rMVPA builds shared train (and, when
present, test) matrices using POSIX shared memory. Task frames carry
voxel or region indices; each worker maps the shared matrices and
extracts the ROI on demand. Model fitting still allocates ordinary
private working memory, and the controller still owns its original rMVPA
objects. Shard therefore reduces data transport; it does not make the
complete analysis zero-copy or eliminate all per-worker memory.

| backend | task_payload | advantages | tradeoffs |
|:---|:---|:---|:---|
| default | controller-extracted ROI matrices | no optional dependency; broad dataset support | serialization and transient ROI payloads |
| shard | ROI indices plus shared-memory handles | small task frames; worker-side extraction; shared source matrix | experimental; setup/cleanup; POSIX platform and supported datasets |
| auto | shard when setup succeeds; otherwise default | convenient opportunistic acceleration | fallback can hide which path actually ran |

rMVPA data-transport backends. {.table}

Use an explicit backend in a production receipt. `auto` is convenient
for an interactive searchlight, but it may fall back to `default` when
shard is absent or the dataset is unsupported. At present shard supports
volumetric, surface, and clustered datasets, but not multibasis image
datasets. High-level runners clean shared segments on normal exit.

``` r

run_with_shared_data <- function(model_spec, radius = 3) {
  if (!requireNamespace("shard", quietly = TRUE)) {
    stop("Install shard before requesting the shared-memory backend.")
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 4)

  run_searchlight(
    model_spec,
    radius = radius,
    backend = "shard"
  )
}
```

Shard is most promising when the source matrices are large, ROIs
overlap, or the same data would otherwise be serialized repeatedly. Its
setup cost can lose on a tiny or single-task analysis. Benchmark both
paths with the same seed, regions, model, workers, and batch size.

## Why can Feature RSA still be slow?

Parallel transport is only one part of
[`feature_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
runtime. For PLS and PCR, the component selector determines how many
numerical fits occur inside every outer cross-validation fold:

| selection | approximate_inner_work | interpretation |
|:---|:---|:---|
| blocked | one fit per training block | block-level held-out error |
| loo | one fit per training observation | observation-level held-out error |
| pve | no validation refits | in-sample explained-variance heuristic |
| max | no validation refits | fixed component count |

Component-selection work inside each Feature RSA outer fold. {.table}

`"loo"` remains the compatibility default. When rows within runs or
sessions are dependent, `ncomp_selection = "blocked"` is usually both
cheaper and a better match to the independence structure. It requires an
observation-aligned `block_var` in
[`feature_rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md)
and at least two training blocks in every outer fold. Inner centering
and scaling are learned without the held-out block. `"pve"` and `"max"`
are cheaper still, but they change the tuning estimand rather than
accelerating the same validation scheme. See
[`vignette("Feature_RSA")`](https://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA.md)
for the selection contract.

rMVPA now streams held-out segment errors, uses compact matrix-level
PLS/PCR fits, and computes the final out-of-fold geometry once. A worker
also drops the original quadratic similarity matrix `S` after its
derived feature targets have been cached, and transports one target
matrix plus fold indices instead of overlapping target slices for every
fold. These changes reduce fitting and serialization overhead for both
the default and shard backends. Shard still matters for the neural
source matrix; it does not eliminate private coefficient, prediction, or
validation work inside each worker.

The source benchmark at `inst/benchmarks/feature_rsa_hotpaths.R`
validates prediction and component-selection parity before reporting
repeated absolute times, compact-fit sizes, and worker-payload sizes.
Treat its timings as a machine-specific characterization, not a
universal threshold.

## How do OpenMP and BLAS interact?

The current rMVPA C source contains no OpenMP parallel regions. That
does not make every rMVPA analysis single-threaded: an optional
classifier, numerical library, BLAS, or another dependency may use
native threads. The active native code path, rather than a package’s
mere link to a runtime, determines whether `OMP_NUM_THREADS` matters.

Set thread limits before R starts. The OpenMP specification defines
`OMP_NUM_THREADS` as an initial control for OpenMP parallel regions;
BLAS implementations have separate controls. Changing a variable after a
runtime has initialized may be too late or runtime-specific.

A useful capacity warning is

``` math
  \text{potential runnable threads}
  \approx \text{Future workers} \times
  \text{maximum native threads per worker}.
```

This is deliberately conservative. OpenMP and BLAS pools need not be
nested or active at the same instant, while some models add another
independent pool. The point is to budget the whole process tree, not to
multiply every setting blindly and call the product a measurement.

Start with one native thread per worker. Increase a native pool only
after a profile shows time in that library, and reduce the Future worker
count so the combined job stays within its CPU and memory allocation. Do
not use a few successful multicore runs with `OMP_NUM_THREADS = 1` as
proof of fork safety: the relevant native runtime may have initialized
earlier, or the test workload may never have entered a parallel region.

## An executable parity check

The following small fixture uses sequential execution as its numerical
reference. It also checks shard when the optional package is installed.
This is a correctness example, not a timing benchmark.

``` r

example <- gen_sample_dataset(
  D = c(4, 4, 4), nobs = 24, nlevels = 2, blocks = 3
)
crossval <- blocked_cross_validation(example$design$block_var)
model_spec <- mvpa_model(
  model = load_model("sda_notune"),
  dataset = example$dataset,
  design = example$design,
  model_type = "classification",
  crossval = crossval,
  return_predictions = FALSE
)

labels <- integer(length(example$dataset$mask))
inside <- which(as.vector(example$dataset$mask) > 0)
labels[inside] <- rep_len(1:4, length(inside))
region_mask <- NeuroVol(labels, space(example$dataset$mask))

reference <- local({
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::sequential)
  run_regional(
    model_spec, region_mask,
    backend = "default", preflight = "off"
  )
})

reference$performance_table
#> # A tibble: 4 × 3
#>   roinum Accuracy     AUC
#>    <int>    <dbl>   <dbl>
#> 1      1    0.5   -0.0556
#> 2      2    0.667  0.222 
#> 3      3    0.458  0.0694
#> 4      4    0.458 -0.0278
```

For parallel validation, keep the fixture fixed and compare the complete
canonical result vector, not just its length, class, or average
accuracy. Run each thread configuration in a fresh process with a
timeout so a low-level hang becomes an observed `timeout` rather than a
stalled notebook.

## What did the benchmark show?

The checked-in benchmark used a 72-observation, `28 x 28 x 28` regional
classification fixture with 18 regions, five batches of at most four
regions, and `sda_notune`. Non-sequential plans configured two workers.
OpenMP and BLAS were fixed at one thread. Each scheduler/data-backend
pairing ran in three fresh R processes from a clean installed rMVPA
0.1.3 tarball.

| Future plan  | Data backend | Median (s) | Range (s)     | Task frames (MiB) |
|:-------------|:-------------|-----------:|:--------------|------------------:|
| sequential   | default      |      5.341 | 5.193–5.523   |            12.444 |
| sequential   | shard        |      4.014 | 3.925–4.108   |             0.101 |
| multisession | default      |      8.577 | 8.459–8.629   |            12.444 |
| multisession | shard        |      6.315 | 6.237–6.480   |             0.101 |
| multicore    | default      |      6.070 | 6.012–6.178   |            12.444 |
| multicore    | shard        |      3.650 | 3.494–4.058   |             0.101 |
| mirai        | default      |     13.737 | 12.690–15.675 |            12.444 |
| mirai        | shard        |     10.891 | 10.022–11.018 |             0.101 |

Local fresh-process characterization on an Apple M3 Max (36 GB). Timings
are descriptive, not thresholds. {.table}

All eight paths produced identical canonical result vectors in every
fixture repetition (`max_abs_result_error = 0`). The default frames
totaled 13,048,800 bytes per run; shard frames totaled 105,712 bytes, a
123-fold reduction for this design. This is direct evidence about the
transport structure measured by
[`object.size()`](https://rdrr.io/r/utils/object.size.html), not a claim
that physical RAM fell by exactly the same factor.

Shard reduced the median time for sequential, multicore, and (slightly)
mirai in this short local run, but it was slower with multisession. The
three-repetition timing ranges overlap for every pairing, and startup
dominates the separate-session plans. Multicore’s local timing does not
override its safety constraints. Mirai’s fresh-process overhead here
does not predict a long-lived mirai deployment. Re-run the driver on the
target system before choosing on speed.

The receipt intentionally omits process-tree RSS as a publication claim.
Summed RSS can double-count shared and copy-on-write pages, while socket
and mirai workers may be detached from the controller’s descendant tree.
Frame payload size, exact output parity, and separately reported timing
phases were the stable measurements.

The reproducibility driver is installed at
`inst/benchmarks/parallel_runtime_benchmark.R`; the more general grid
and its field definitions are in `scripts/sweep_parallel_runtime_grid.R`
and `scripts/README_parallel_runtime_grid.md` in the source repository.

## Practical decisions

Use this order of operations:

1.  Run the complete analysis sequentially on a small fixture and save
    its canonical result.
2.  Determine the CPU and memory allocation with
    [`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
    and the scheduler, not the host’s advertised core count.
3.  Select `multisession`, start with native thread pools at one, and
    verify that the intended number of batches and workers actually ran.
4.  Compare `default` and `shard` explicitly. Keep shard only if its
    smaller task frames or elapsed time matter for the real workload.
5.  For PLS/PCR Feature RSA, choose the component-selection estimand
    before tuning workers. Prefer `"blocked"` when the block is the
    independence unit; exact `"loo"` requires one inner fit per training
    observation.
6.  Consider mirai for its scheduler and lifecycle properties, not as a
    synonym for shared memory.
7.  Increase one axis at a time. Recheck complete numerical results
    after every change.
8.  Use multicore only when the full native stack is known to be
    fork-safe; do not combine it casually with OpenMP or another
    multithreaded runtime.

The underlying runtime references are the [Future plan
documentation](https://future.futureverse.org/reference/plan.html),
[Future multicore
warning](https://future.futureverse.org/reference/multicore.html),
[future.mirai documentation](https://future.mirai.futureverse.org/),
[furrr gotchas](https://furrr.futureverse.org/articles/gotchas.html),
[parallelly core
detection](https://parallelly.futureverse.org/reference/availableCores.html),
[R Core fork
guidance](https://stat.ethz.ch/R-manual/R-devel/library/parallel/help/mcfork.html),
and the [OpenMP environment-variable
specification](https://www.openmp.org/spec-html/5.0/openmpch6.html).
