# HR‑FReM (Hybrid Randomized Multi‑Scale FReM) — Implementation Plan for rMVPA

Goal: ensemble decoder with multi‑scale randomized parcellations, cluster‑level glmnet fits, voxel back‑projection, and stability scoring; CPU‑friendly, spatially smooth, reproducible.

## Milestones
- **M1: Rcpp kernels** — implement `reduce_features_cpp`, `expand_weights_cpp` (counts-aware, optional OpenMP); add unit tests.
- **M2: Parcellation library (R)** — build randomized multi‑scale parcellations; cache `clusters` and `X_red` per entry; stub clustering first, upgrade to `rena_cpp`.
- **M3: HR‑FReM model spec** — constructor + defaults; S3 methods to integrate with `run_searchlight` / `run_regional`.
- **M4: Ensemble training loop** — bootstrap, parcellation draw, glmnet fit (global λ or OOB), back‑project, aggregate weights/stability.
- **M5: Prediction helper** — apply final voxel weights to new data; basic checks.
- **M6: Docs & examples** — roxygen for new functions; short vignette snippet; note CPU/memory tips.

## Module Sketch
- `src/hr_frem_ops.cpp`: RcppArmadillo kernels (reduce/expand); later add `rena_cpp`.
- `R/hr_frem_library.R`: `build_parcellation_library(X, K_levels, M_per_scale, rena_fun, ...)` with cached `X_red`.
- `R/hr_frem_model.R`: constructor (`hr_frem_model`), args (K_levels, M_per_scale, B, alpha, lambda_mode ∈ {global,oob}, lambda_star, bag_weighting, seed, library=NULL, return_predictions=FALSE).
- `process_roi.hr_frem_model`: phases:
  - (optional) global λ* via `cv.glmnet` on one `X_red`.
  - loop over bags: stratified bootstrap, sample parcellation, slice `X_red`, glmnet fit, OOB λ selection if needed, back‑project weights, accumulate `w_accum`, `stability_cnt`, OOB acc.
  - aggregate: `w_final = w_accum/weight_sum`; `stability = stability_cnt/B`; return classification_result + diagnostics.
- `R/hr_frem_predict.R`: predict with voxel weights (`type = response|class`).

## Testing Checklist
- Rcpp correctness: cluster means vs R aggregate; expand(reduce(X)) preserves voxel sums; NA/degenerate cluster guards.
- Synthetic data: planted sparse blobs → AUC > chance; stability peaks at signal voxels.
- Path coverage: lambda_mode = "global" and "oob"; bag_weighting on/off; determinism under fixed seed.
- Memory/size sanity on moderate p (e.g., p≈20k) to confirm speed.

## Optional Enhancements (later)
- Replace placeholder clustering with spatial `rena_cpp` (ReNA/Ward‑style, randomized, adjacency-aware).
- Disk‑backed library (fst/qs) to avoid keeping all `X_red` in RAM.
- Parallel bags via `future` with proper RNG streams.

## Integration Notes
- Update `DESCRIPTION` (Imports: Rcpp, glmnet; LinkingTo: Rcpp, RcppArmadillo); `NAMESPACE` useDynLib registration; run `Rcpp::compileAttributes()`.
- Keep default `return_predictions=FALSE` to avoid heavy payloads; store per‑bag diagnostics only if requested.
