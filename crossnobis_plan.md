# Implementation Plan – Crossnobis Support for MS-ReVE / Contrast-RSA

## 0 Paper Summary & Motivation
The cross-nobis (a.k.a. *cross-validated Mahalanobis*) distance provides an **unbiased** estimate of the true squared pattern distance by **multiplying pattern differences that were estimated from *independent* data partitions**.  Ordinary (same-fold) distances are positively biased because noise is multiplied by itself.

Formal definition for condition pair \(k=(i,j)\) with \(M\) partitions (Eq. 5, Diedrichsen & Kriegeskorte 2017):
\[
  \tilde d_k 
    = \frac{1}{M(M-1)\,P}\sum_{m\neq n} \hat\delta_{k,m}^{\mathsf T}\,\hat\delta_{k,n},
  \qquad \hat\delta_{k,m}=\hat b_{i,m}-\hat b_{j,m},
\]
where \(P\) is the number of voxels/channels.  Equivalent "single-fold-vs-rest" algebra can also be used.

Key properties
* **Unbiased**: \(\mathbb E[\tilde d_k]=\delta_k^{\mathsf T}\delta_k/P\).
* Can be **negative** when the true distance is ~0 – keep these values, do not truncate.
* Variance is \(\tfrac{M}{M-1}\) times larger than the biased estimator but decreases with more partitions.
* **Whitening**: pre-multiply patterns by \(\mathbf W=\Sigma_{\text{noise}}^{-1/2}\) to obtain Mahalanobis distances; crucial for fMRI.

---

(Sections 1–9 below incorporate this theory.)

## 1 High-Level Requirements
*Produce unbiased squared distances …*  
*(unchanged from previous draft)*

// ... existing content of sections 1–9 from crossnobis.md retained unchanged ...

## 9 Timeline (suggestion)
| Week | Deliverable |
|------|-------------|
| 1 | Helper function + unit tests pass |
| 2 | Integration into train_model + searchlight plumbing |
| 3 | Documentation, vignette, example rebuild |
| 4 | Profiling & optional Rcpp optimisation |

---

## Appendix A  Detailed Implementation Table & Sanity-Check Code

### A.1 Step-wise Implementation Table (adapted from user note)
| Step | What to compute | Suggested R objects / code changes |
|------|-----------------|------------------------------------|
| **A** | *Return fold-wise condition means.*  | `compute_crossvalidated_means_sl(return_folds = TRUE)` returns a list with `mean = U_hat_sl`, `folds = U_folds` (K × V × M). |
| **B** | *Pre-whiten each fold.* | Within the same helper: if `W` provided, `U_folds[,,m] %*% W`.
| **C** | *Cross-fold inner products only for m ≠ n.* | New helper `crossnobis_distance(U_folds, P)` computes distance vector; pseudo-code shown below. |
| **D** | *Expose to RSA regression.* | `train_model.contrast_rsa_model()` uses this distance vector when `estimation_method == "crossnobis"`.  `include_vec` still applied. |
| **E** | *Allow negative outputs.* | No truncation; document in Roxygen. |
| **F** | *Unit tests.* | Noise-only simulation: mean distance ≈ 0; biased Euclidean positive. |
| **G** | *Performance tips.* | Vectorise across pairs; optionally migrate to Rcpp. |

```r
crossnobis_distance <- function(U_folds) {
  K <- dim(U_folds)[1]; V <- dim(U_folds)[2]; M <- dim(U_folds)[3]
  if (M < 2) rlang::abort("Need at least two partitions for cross-nobis.")
  pair_mat <- utils::combn(K, 2)
  D <- ncol(pair_mat)
  out <- numeric(D)
  P <- V
  for (p in seq_len(D)) {
    i <- pair_mat[1, p]; j <- pair_mat[2, p]
    deltas <- U_folds[i, , ] - U_folds[j, , ]   # V × M
    ip <- tcrossprod(t(deltas))                 # M × M inner products
    diag(ip) <- 0
    out[p] <- sum(ip) / (P * M * (M - 1))
  }
  names(out) <- paste0(dimnames(U_folds)[[1]][pair_mat[1, ]], "_vs_", 
                       dimnames(U_folds)[[1]][pair_mat[2, ]])
  out
}
```

### A.2 Quick Sanity-Check
```r
set.seed(1)
K <- 4; V <- 30; M <- 8
true <- matrix(rnorm(K*V), K, V)          # ground-truth patterns
noise <- array(rnorm(K*V*M, sd = 3), c(K, V, M))
U_folds <- true + noise                   # simulated noisy means per fold
cross_d <- crossnobis_distance(U_folds)
mean(cross_d)   # ≈ 0 for pure noise case
```

## 10 Checkable Implementation Tickets

This section breaks down the implementation plan into granular, checkable tasks for the specified R files and the new helper function.

Always check off each task as you complete it.

### New Helper Function (e.g., in `R/crossnobis_helpers.R` or similar)

*   [ ] **Implement `compute_crossnobis_distances_sl` function:**
    *   [ ] Define function signature: `compute_crossnobis_distances_sl(U_folds, P_voxels)`.
        *   `U_folds`: K x V x M array of condition means per fold (K=conditions, V=voxels, M=folds).
        *   `P_voxels`: Number of voxels (V).
    *   [ ] Handle edge case: If `M < 2`, abort with an informative message.
    *   [ ] Generate all unique unordered condition pairs (e.g., using `utils::combn(K, 2)`).
    *   [ ] For each condition pair `k = (i,j)`:
        *   [ ] Initialize `sum_cross_prods_k = 0`.
        *   [ ] For each pair of distinct folds `(m, n)` where `m != n`:
            *   [ ] Calculate `delta_k_m = U_folds[i,,m] - U_folds[j,,m]`.
            *   [ ] Calculate `delta_k_n = U_folds[i,,n] - U_folds[j,,n]`.
            *   [ ] Add `crossprod(delta_k_m, delta_k_n)` to `sum_cross_prods_k`.
        *   [ ] Calculate `d_crossnobis_k = sum_cross_prods_k / (P_voxels * M * (M - 1))`.
    *   [ ] Return a named numeric vector of these `d_crossnobis_k` values, with names like "condA_vs_condB", ordered to match lower-triangle vectorization.
    *   [ ] Handle potential NAs in `U_folds` (e.g., if a condition is missing in a fold): ensure resulting distances involving such conditions are NA, or relevant pairs are skipped.
    *   [ ] Add Roxygen documentation for this new function.

### `R/compute_cv_means.R` (Changes to `compute_crossvalidated_means_sl`)

*   [ ] **Argument and Return Value Modification:**
    *   [ ] Add a new boolean argument `return_folds = FALSE` to `compute_crossvalidated_means_sl`.
    *   [ ] If `return_folds = TRUE`, change the return value to a list: `list(mean_estimate = U_hat_sl, fold_estimates = U_folds_array)`, where `U_folds_array` is (K x V x M).
    *   [ ] If `return_folds = FALSE` (default), return `U_hat_sl` as before (or `mean_estimate` from the list).
*   [ ] **Per-Fold Mean Calculation (when `return_folds = TRUE`):**
    *   [ ] Inside the loop over folds, store the `fold_means_mat_processed` (K_fold x V matrix for conditions present in the current fold) before re-indexing to `full_fold_means`.
    *   [ ] Collect these `fold_means_mat_processed` (or their `full_fold_means` equivalent containing NAs for absent conditions) into the `U_folds_array` (K x V x M).
*   [ ] **Whitening Logic (when `return_folds = TRUE` and `whitening_matrix_W` is provided):**
    *   [ ] Ensure that if `estimation_method == "crossnobis"` (as passed to `compute_cv_means`) and `whitening_matrix_W` is provided, the whitening (`fold_means_mat %*% whitening_matrix_W`) is applied to the patterns that will form `U_folds_array`. (Current logic already does this for `fold_means_mat_processed`).
*   [ ] **Documentation:**
    *   [ ] Update Roxygen for `compute_crossvalidated_means_sl` detailing `return_folds` and the new list return type.

### `R/contrast_rsa_model.R`

*   **`contrast_rsa_model` Constructor:**
    *   [ ] Re-enable `estimation_method = "crossnobis"` in `match.arg()`. (Remove the `rlang::abort` for it).
    *   [ ] Add argument `whitening_matrix_W = NULL` to the function signature.
    *   [ ] Store `whitening_matrix_W` in the model specification object (e.g., `x$whitening_matrix_W <- whitening_matrix_W`).
    *   [ ] Update Roxygen:
        *   [ ] Reflect "crossnobis" as a valid `estimation_method`.
        *   [ ] Document the new `whitening_matrix_W` argument, its purpose for crossnobis, and that it's passed to `compute_crossvalidated_means_sl`.

*   **`train_model.contrast_rsa_model` Method:**
    *   [ ] **Initialization:**
        *   [ ] Declare `U_hat_for_delta_calc` (will hold the K x V matrix for Δ projections).
        *   [ ] Declare `dvec_sl` (will hold the vector of (dis)similarities for regression).
    *   [ ] **Conditional Logic for `estimation_method`:**
        *   [ ] **If `obj$estimation_method == "crossnobis"`:**
            *   [ ] Call `compute_crossvalidated_means_sl` with:
                *   `estimation_method = "crossnobis"` (this tells `compute_cv_means` to handle whitening internally if `W` is passed).
                *   `whitening_matrix_W = obj$whitening_matrix_W`.
                *   `return_folds = TRUE`.
                *   Store result in `cv_outputs`.
            *   [ ] Assign `U_hat_for_delta_calc <- cv_outputs$mean_estimate`.
            *   [ ] Let `U_folds_data <- cv_outputs$fold_estimates`.
            *   [ ] `P_voxels <- ncol(sl_data)`.
            *   [ ] `dvec_sl <- compute_crossnobis_distances_sl(U_folds_data, P_voxels)`.
        *   [ ] **Else (for "average", "L2_norm"):**
            *   [ ] `U_hat_for_delta_calc <- compute_crossvalidated_means_sl(sl_data, mvpa_des, cv_spec, obj$estimation_method, whitening_matrix_W = NULL)` (Note: `compute_cv_means` correctly ignores `W` if its own `estimation_method` is not "crossnobis").
            *   [ ] `G_hat_sl <- U_hat_for_delta_calc %*% t(U_hat_for_delta_calc)`.
            *   [ ] `dvec_sl <- G_hat_sl[lower.tri(G_hat_sl)]`.
    *   [ ] **Common Post-Processing for `dvec_sl` (after if/else block):**
        *   [ ] If `obj$estimation_method != "crossnobis"`, apply `include_vec` logic to `dvec_sl` based on `G_hat_sl` as currently done.
        *   [ ] If `obj$estimation_method == "crossnobis"`, apply `include_vec` logic directly to the `dvec_sl` output by `compute_crossnobis_distances_sl` (ensure mapping from KxK `is_intra_block_pair` to K(K-1)/2 `dvec_sl` elements is correct).
    *   [ ] **Voxel Projections (Δ_sl):**
        *   [ ] Calculate `Delta_sl <- t(U_hat_for_delta_calc) %*% C_ord` using the consistently derived `U_hat_for_delta_calc`.
    *   [ ] (Rest of the RSA regression, metric calculation proceeds using `dvec_sl` and `beta_sl`, and `Delta_sl` as appropriate).
*   **`run_searchlight.contrast_rsa_model` Method:**
    *   [ ] No specific changes required here if `whitening_matrix_W` is stored in and accessed from `model_spec` by `train_model`.

### `R/msreve_design.R`

*   [ ] No changes anticipated for this file based on the current `crossnobis_plan.md`.
