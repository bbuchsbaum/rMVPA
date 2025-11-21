# kernel_searchlight: a pluggable searchlight engine

## Motivation
- Current `run_searchlight()` mixes three concerns: (1) where to sample (grid/random), (2) which voxels belong to each ROI (shape), and (3) what computation to run (MVPA classifier, RSA, REMAP-RRR). That tight coupling makes it hard to add anisotropic shapes, importance sampling, or swap in different metrics without duplicating code paths.
- A light abstraction layer lets us re-use the same sampling/geometry machinery across classic decoding, RSA, crossnobis, and REMAP-RRR, while preserving the existing API for backward compatibility.

## Core idea
Expose four small, composable hooks with safe defaults:
1) **sampler**: yields ROI centers (grid, Poisson-disk, random mask, surface nodes).
2) **shape_bank**: defines voxel offsets (spheres, ellipsoids, GM-aligned variants) plus optional soft weights.
3) **materializer**: turns a `(center, shape_id)` pair into concrete voxel indices/weights (handles rotations, boundary clipping, surface/volume differences).
4) **kernel**: the computation on the extracted ROI data (e.g., SVM/LDA decoding, RSA correlations, crossnobis, REMAP-RRR’s `process_roi`). It may return one or many metrics (vector/list) plus optional diagnostics; aggregation is handled separately so multi-metric kernels remain supported.

`kernel_searchlight()` stitches these together; `run_searchlight()` can simply delegate to it with the current grid + spherical + SVM defaults, so existing scripts keep working.

## How it improves abstraction and flexibility
- **Pluggable computation**: Kernels become first-class objects; adding RSA, crossnobis, or REMAP-RRR is just picking a different kernel, not rewriting the searchlight loop.
- **Geometry decoupled from stats**: Shapes and materializers handle anisotropy/GM alignment; kernels stay agnostic to voxel layout.
- **Sampling strategies swap cleanly**: Even coverage (Poisson-disk), random subsets for speed, or surface-based samplers can be dropped in without touching kernels.
- **Backward compatible**: Defaults mirror current behavior; legacy calls to `run_searchlight()` are unaffected.
- **Future-proof**: New ideas—adaptive refinement, weighted aggregators, GPU-friendly soft shapes—fit by adding hook implementations, not new monoliths.

### On "metrics" as a separate concern
- Metrics often sit inside kernels (e.g., RSA may emit multiple similarity stats; REMAP-RRR returns accuracy, AUC, adaptor diagnostics). Pulling "metric" out as a universal hook is attractive, but many kernels bundle task-specific logic (e.g., class labels, crossvalidation folds) that a generic metric stage would need to re-learn.
- Pragmatic compromise: kernels return a named vector/list of metrics; the **aggregator** chooses which metric(s) to project back to voxels/ROIs (e.g., pick `accuracy`, `auc`, or a user-specified `metric_name`). This keeps multi-metric support without forcing every kernel into a common metric interface.

**Standardized kernel return format**: To prevent aggregator complexity from handling arbitrary structures, enforce a simple contract:
```r
kernel_result <- list(
  metrics = c(accuracy = 0.75, auc = 0.82, ...),  # named numeric vector (required)
  diagnostics = list(...)                         # optional, arbitrary structure
)
```
This ensures aggregators can always assume `result$metrics` exists as a flat numeric vector, while diagnostics (like REMAP's `diag_by_fold`) live in a separate namespace. Multi-metric kernels work uniformly without complex introspection.
## Parallel plan for regional methods (mvpa_regional & friends)
- Introduce a sibling entry point: `kernel_regional()` that mirrors the hook pattern but replaces sampler/shape with a **roi_provider** (atlas labels, mask list, parcels, surfaces).
- Hooks: `roi_provider` → yields ROI index sets; `roi_materializer` → optional weighting/cleaning (exclude NaNs, GM-only); `kernel` → same kernels as searchlight (SVM/LDA, RSA, crossnobis, REMAP-RRR); `aggregator` → how to combine multi-metric outputs per ROI.
- Rationale: shared kernels keep stats code single-sourced; switching between parcel-based and searchlight analyses becomes swapping providers, not rewriting pipelines.
- Nice side effects: cache per-ROI covariance/whitening across kernels, balance workload by ROI size, enable batched/parallel execution identical to searchlight backends.

## Design considerations and refinements

### Shape bank vs materializer: simplification option
Currently split into two hooks for maximum flexibility:
- **shape_bank**: defines abstract voxel offsets
- **materializer**: instantiates concrete indices/weights from (center, shape)

**Alternative**: Merge into single `shape_provider(center, params) → (indices, weights)`
- **Pro**: Simpler mental model; one function call instead of two
- **Con**: Loses ability to pre-define/cache shapes independently

**Recommendation**: Keep the separation but provide convenience factories for common cases:
```r
spherical_provider(radius = 8)           # bundles shape + materializer
ellipsoidal_provider(radii = c(8,8,4))   # anisotropic
surface_geodesic_provider(radius_mm = 15)
```
This gives power users full control while simplifying the typical workflow.

### Weighting support: incremental adoption
Materializers can return soft weights (GM probability, distance-based falloff), but not all kernels support weighted observations:
- **SVM/glmnet**: Already support `weights` parameter
- **RSA**: Would need weighted correlation (doable but not implemented)
- **REMAP-RRR**: Currently unweighted; would require changes to joint whitening

**Approach**:
1. Start with binary (0/1) weights for inclusion/exclusion
2. Add `supports_weights` metadata flag to kernel specs
3. Soft weighting becomes incremental enhancement as kernels gain capability
4. Aggregator can warn if soft weights provided to non-supporting kernel

### Surface data: first-class support
The materializer abstraction makes native surface support cleaner than adapters:
- `surface_geodesic_materializer(center_vertex, radius_mm)` → vertex indices/weights
- Kernels receive vertex data identically to voxel data
- Since rMVPA already integrates `neurosurf`, native surface materializers are worthwhile

### Execution strategy: 5th hook for parallelization
Different kernels have different resource profiles; consider adding an optional `executor` hook:
```r
kernel_searchlight(
  sampler = grid_sampler(...),
  materializer = spherical(...),
  kernel = remap_rrr_kernel(...),
  executor = future_executor(batch_size = 100)  # controls parallelization
)
```
**Benefits**:
- Light kernels (SVM) use fine-grained parallelism
- Heavy kernels (REMAP-RRR with `return_adapter=TRUE`) use coarse batching
- Surface searchlights can choose node-parallel vs batch-parallel
- Future GPU kernels get appropriate scheduling

**Default**: Current `mvpa_iterate` behavior for backward compatibility

### Aggregator specification
The aggregator contract needs explicit definition to handle multi-metric outputs, errors, and combining strategies:

```r
aggregator(
  results,        # tibble: roi_id, result, indices, error, error_message, ...
  metric_names,   # which metrics to extract from kernel results (default: all)
  combine_fn,     # for randomized/resampled: mean, median, pool_randomized
  dataset         # for building spatial maps (NeuroVol/NeuroSurface)
) → list(
  results = list(accuracy = NeuroVol(...), auc = NeuroVol(...)),
  n_voxels = ...,
  active_voxels = ...,
  metrics = c("accuracy", "auc")
)
```

**Error handling**: When `error=TRUE` for some ROIs, aggregator:
1. Warns about failed ROI count
2. Assigns `NA` to those voxels in output maps
3. Reports `active_voxels` (successful) separately from `n_voxels` (total)

Current `combine_randomized`, `combine_standard`, `pool_randomized` become aggregator implementations.

### Degenerate cases and error policies
**Empty ROIs** (materializer returns 0 voxels):
- Kernel receives empty data, should return `list(metrics = NA_real_, diagnostics = NULL)`
- Aggregator propagates `NA` to spatial map for that location

**Insufficient data** (e.g., REMAP-RRR `<2 paired items`):
- Current: hard error with `error=TRUE` in result tibble
- Alternative: warning + graceful `NA` return
- **Decision**: Kernel-specific policy; provide both `stop()` and `NA`-return patterns
  - Critical failures (missing test data) → error
  - Degenerate cases (too few observations) → NA with warning

**Memory limits**:
- Large searchlights with `return_adapter=TRUE` or heavy diagnostics can OOM
- **Recommendation**: Document memory requirements; consider `drop_diagnostics` parameter for searchlight mode
- Diagnostics should be opt-in for searchlight, default for regional

## Implementation roadmap

### Phase 1: Core infrastructure (backward compatibility focus)
**Goal**: Establish hook interfaces without breaking existing code

1. Define hook interfaces:
   - `sampler(dataset, ...) → iterator of centers`
   - `materializer(center, dataset, ...) → (indices, weights)`
   - `kernel(roi_data, model_spec, ...) → list(metrics = c(...), diagnostics = list(...))`
   - `aggregator(results, metric_names, combine_fn, dataset) → searchlight_result`

2. Implement default hooks mirroring current behavior:
   - `grid_sampler()` → current grid logic
   - `spherical_materializer(radius)` → current spherical ROI extraction
   - `mvpa_kernel()` → wraps current MVPA model fitting
   - `standard_aggregator()` → current result combination

3. Create `kernel_searchlight()` as orchestrator:
   ```r
   kernel_searchlight(dataset, model_spec,
                     sampler = grid_sampler(),
                     materializer = spherical_materializer(radius = 8),
                     kernel = mvpa_kernel(),
                     aggregator = standard_aggregator())
   ```

4. Refactor `run_searchlight()` to delegate:
   ```r
   run_searchlight.mvpa_model <- function(model_spec, radius = 8, ...) {
     kernel_searchlight(model_spec$dataset, model_spec,
                       materializer = spherical_materializer(radius),
                       ...)
   }
   ```

5. **Validation**: Run full test suite; all existing tests must pass unchanged

### Phase 2: Kernel extraction (code consolidation)
**Goal**: Eliminate duplication by extracting model-specific logic into kernels

1. Extract `mvpa_kernel()` from current MVPA model fitting loop
2. Extract `rsa_kernel()` from RSA searchlight implementation
3. Extract `remap_rrr_kernel()` from `process_roi.remap_rrr_model()`
4. Extract `crossnobis_kernel()` if applicable
5. Each kernel should handle its own cross-validation, metric computation, diagnostics
6. **Validation**: Feature parity checks; compare old vs new results numerically

### Phase 3: New capabilities (extensibility)
**Goal**: Add functionality impossible or awkward in current design

1. **Sampling strategies**:
   - `poisson_disk_sampler(min_distance)` → even coverage
   - `random_sampler(n_samples)` → subset for speed
   - `surface_node_sampler()` → for surface data

2. **Shape variations**:
   - `ellipsoidal_materializer(radii = c(x, y, z))` → anisotropic
   - `gm_weighted_materializer(radius, gm_prob_vol)` → GM-based soft weights
   - `surface_geodesic_materializer(radius_mm)` → surface patches

3. **Convenience factories**:
   - `spherical_provider(radius)` → bundles shape + materializer
   - `searchlight_preset(type = "fast")` → pre-configured pipelines

4. **Documentation**: Vignette on "Custom Kernels and Materializers"

### Phase 4: Regional unification
**Goal**: Share kernels between searchlight and regional analyses

1. Implement `kernel_regional()` with parallel structure:
   ```r
   kernel_regional(dataset, model_spec,
                  roi_provider = atlas_provider(mask),
                  roi_materializer = identity_materializer(),
                  kernel = same_kernels_as_searchlight,
                  aggregator = regional_aggregator())
   ```

2. Refactor `run_regional()` to delegate to `kernel_regional()`

3. **Benefits realized**:
   - Single kernel codebase for both analysis modes
   - Easy A/B testing: same model, searchlight vs parcels
   - Cached preprocessing (covariance, whitening) across kernels

### Phase 5: Advanced features (optimization)
**Goal**: Performance and GPU preparation

1. Add `executor` hook for custom parallelization strategies
2. Batch-aware kernels for memory efficiency
3. Caching layer for expensive preprocessing (GM masks, covariance matrices)
4. Profiling-guided defaults (auto-detect optimal batch size)
5. GPU-ready kernel interface (data layout, memory pinning)

## Testing strategy

### Regression tests
- **Numeric equivalence**: New kernel-based code must produce bit-identical results to current implementation
- **API compatibility**: All existing `run_searchlight()` calls work without modification
- **Edge cases**: Empty ROIs, degenerate data, single-voxel masks

### New functionality tests
- Custom samplers/materializers produce expected ROI counts and geometries
- Multi-metric kernels correctly populate all output maps
- Error handling: graceful degradation when kernels fail on subset of ROIs
- Memory: large searchlights with diagnostics don't OOM (or fail predictably)

### Integration tests
- Searchlight → regional workflow with same kernel
- Surface data end-to-end
- Parallel execution across backends (sequential, future, foreach)
