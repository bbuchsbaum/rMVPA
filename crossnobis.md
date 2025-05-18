# Implementation Plan – Crossnobis Support for MS-ReVE / Contrast-RSA

**Status (May 2025)**  
`estimation_method = "crossnobis"` is currently disabled in `contrast_rsa_model()`.  The code only whitened fold-means and did **not** compute the unbiased squared Euclidean distances required by Eq. (5) of Diedrichsen & Kriegeskorte (2017).  This document specifies the steps to implement Crossnobis correctly, end-to-end.

---
## 1 High-Level Requirements
* Produce unbiased squared distances \(\tilde d_k\) for every condition pair \(k = (i,j)\) inside each searchlight / ROI.
* Follow Eq. (5):
  \[ \tilde d_k = \frac{1}{M(M-1) P} \sum_{m=1}^M \sum_{\substack{n=1\\ n\ne m}}^M \hat\delta_{k,m}^{\top}\, \hat\delta_{k,n} \]
  where \(M\) = # partitions, \(P\)=# voxels, \(\hat\delta_{k,m} = \hat\mu_{i,m}-\hat\mu_{j,m}\).
* Optional noise-whitening via a pre-computed matrix \(\mathbf W = \Sigma_{noise}^{-\tfrac12}\).
* Integrate seamlessly with current RSA regression pipeline, preserving β/Δ logic and output metrics.

---
## 2 API Surface Changes
| Component | Change |
|-----------|--------|
| `contrast_rsa_model()` | Re-allow `estimation_method = "crossnobis"`; add argument `whitening_matrix_W = NULL` that is stored in the model spec. |
| `compute_crossvalidated_means_sl()` | **No change** – new function will be used. |
| **New helper** `compute_crossnobis_distances_sl()` | Returns a *named* vector `d_crossnobis` of length \(K(K-1)/2\). |
| `train_model.contrast_rsa_model()` | Branch: if `obj$estimation_method == "crossnobis"`, call the new helper to obtain `dvec_sl` directly **(skip Ĝ, lower-tri step)**.  Still compute `U_hat_sl` (via existing function) so that Δ and Σ_q β_q Δ_q,v remain available. |
| Searchlight plumbing | Ensure `whitening_matrix_W` is forwarded (either via `...` at `run_searchlight()` or by storing it in the model spec and passing down). |

---
## 3 Algorithm – `compute_crossnobis_distances_sl()`
1. **Inputs**
   * `sl_data` (N × P)
   * `mvpa_design` (with `Y` conditions & `block_var`/cv grouping)
   * `cv_spec` (provides `get_nfolds()` & `train_indices()`)
   * `whitening_matrix_W` (optional)
2. **Preparation**
   * Let `cond_levels` = levels(mvpa_design$Y).
   * Determine `M` folds and allocate accumulator `cross_sum` (length = Kpairs) initialised to 0.
3. **Per-fold means**
   * For each fold `m`:
     * `μ̂_{c,m}` = mean pattern of condition `c` computed **only on training samples of that fold** (same logic as existing mean function).
     * If whitening requested: `μ̂_{c,m} ← μ̂_{c,m} %*% W`.
   * Store the means in an array `M_folds[m, c, p]`.
4. **Cross-products**
   * For every unordered condition pair `(i,j)` (vector index `k`):
     * For all ordered fold pairs `(m,n)` with `m ≠ n`:
       * `δ̂_{k,m} ← μ̂_{i,m} − μ̂_{j,m}`
       * `δ̂_{k,n} ← μ̂_{i,n} − μ̂_{j,n}`
       * `cross_sum[k] += crossprod(δ̂_{k,m}, δ̂_{k,n})`
   * Complexity: O(M² Kpairs P).  For typical `M≤10`, `K≤12`, `P≤500`, feasible.
   * Memory: can compute `δ̂` on the fly; no need to store full array if we loop cleverly.
5. **Normalisation**
   * `d_crossnobis[k] ← cross_sum[k] / (M*(M-1)*P)`
6. **Return**
   * Named vector with names like `"condA_vs_condB"` in lower-tri order to match how contrast RDMs are vectorised elsewhere.

---
## 4 Integration into RSA Regression
* Replace earlier logic building `dvec_sl` from `Ĝ_sl`.
* `include_vec` handling (exclude intra-run pairs) still applies – we set forbidden entries to `NA` *before* regression.
* Predictor matrix `Xmat` (contrast RDMs) is unchanged; we only swap the dependent variable.

---
## 5 Edge-Case Handling
* Missing conditions in some folds → if any `μ̂_{c,m}` undefined, exclude all pairwise distances involving that condition from regression; flag with `na_reason`.
* `M < 2` (e.g., leave-one-out CV) → abort: cannot form cross-products with `m≠n`.
* Near-singular `W` → warn & skip whitening.

---
## 6 Tests & Validation
1. **Unit tests (tests/testthat)**
   * Simulated data with known true pattern differences; verify bias of naïve vs. crossnobis distance.
   * Check invariance to permutation of folds.
   * Confirm equal results when `W = I` vs. no whitening.
2. **Integration test**
   * Run a tiny searchlight with `contrasts_rsa_model(estimation_method="crossnobis")`; ensure no error and reasonable output.

---
## 7 Documentation Updates
* Update Roxygen in all touched functions.
* Add a vignette section explaining the mathematics and when to use crossnobis.

---
## 8 Migration Steps
1. **Phase 1** – Backend logic
   * Implement `compute_crossnobis_distances_sl()` + tests.
   * Wire into `train_model.contrast_rsa_model()`.
2. **Phase 2** – API & examples
   * Re-enable constructor option, add parameter `whitening_matrix_W`.
   * Update searchlight examples.
3. **Phase 3** – Performance optimisation (optional)
   * C++/Rcpp implementation for heavy workloads.

---
## 9 Timeline (suggestion)
| Week | Deliverable |
|------|-------------|
| 1 | Helper function + unit tests pass |
| 2 | Integration into train_model + searchlight plumbing |
| 3 | Documentation, vignette, example rebuild |
| 4 | Profiling & optional Rcpp optimisation |

---
