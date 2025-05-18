# Contrast RSA / MS-ReVE Phase 2 Implementation Plan

This document outlines the features and refinements to be added to the `contrast_rsa_model` functionality based on the original proposal (`signed_rsa_proposal.md`) and the current implementation state.

Keep things modular and simple!
No God functions!
No functions more than 100 lines!

** always check off checkboxes when done **


## Core Engine Enhancements (Section C of Proposal)

- [x] **Advanced Û Estimation (C.1):** Implemented `"L2_norm"` and `"crossnobis"` (requiring user-supplied whitening matrix `W`) options in `compute_crossvalidated_means_sl`. (Further methods like internal Σ estimation for crossnobis from betas could be future enhancements if needed, but current implementation relies on pre-computed `W` derived from GLM residuals as per best practice).
- [x] **RSA Regularization (C.3):** Implement regularization options (`"ridge_hkb"` added for Hoerl-Kennard-Baldwin lambda). (Further options like GCV, fixed lambda, or elastic-net still pending for future enhancement).
- [x] **RSA Collinearity Check (C.3):** Implement the `check_collinearity=TRUE` functionality in the `contrast_rsa_model` constructor to check the rank of the contrast RDM predictor matrix (`Xmat`).
- [x] **Normalized Δ_q Contributions (C.4):** Add an option (`normalize_delta` parameter in `contrast_rsa_model`) to calculate and potentially use normalized voxel contributions (`~Δ_q = Δ_q / ||Δ_q||`).

## Voxel-Level Refinements (Section D of Proposal)

- [ ] **Voxel Reliability Weighting (ρ_q,v) (D.1):** Implement the calculation of voxel contribution reliability across cross-validation folds within `train_model.contrast_rsa_model`. This requires storing fold-specific `U_hat_sl` or `Delta_sl`.

    > **Context: Why Voxel Reliability (ρ_q,v)?**
    >
    > The `β · Δ` map answers: *What direction and how strongly does voxel v contribute to contrast q on average within this searchlight?*
    >
    > Reliability (ρ_q,v) answers the follow-up: *Is that contribution repeatable across different subsets of the data (folds)?*
    >
    > **Intuition:** If Δ_q,v jumps wildly from fold to fold, it's likely driven by noise. ρ scales from 0 (unreliable) to 1 (perfectly stable).
    >
    > **Formula:**
    > \[ \rho_{q,v} = 1 - \frac{ \text{Var}_{\text{folds}}( \Delta_{q,v}^{(\text{fold})} ) }{ \text{Var}_{\text{folds}}( \Delta_{q,v}^{(\text{fold})} ) + \sigma^2_{\text{noise},q,v} } \quad \in [0, 1] \]
    >
    > - Numerator: Empirical fold-to-fold variance of the voxel's contribution.
    > - Denominator: Adds expected noise variance (σ²_noise) to prevent over-penalizing stable signals measured with limited data.
    >
    > **Noise Variance Estimation (σ²_noise,q,v):**
    > This represents the expected variance of Δ_q,v if the true signal were zero. It's derived from within-condition residual variances (σ²_k,v(s)) estimated from the noise whitening/pre-processing stage, weighted by the contrast weights (w_qk) and the number of trials (n_k(s)) per condition k in fold s:
    > \[ \sigma^2_{\text{noise},q,v} = \frac{1}{S} \sum_{s=1}^{S} \sum_{k=1}^{K} \frac{ w_{qk}^2 }{ n_k^{(s)} } \sigma_{k,v}^{2(s)} \]
    > (This requires access to within-condition residual variances, which might need to be passed through or estimated differently depending on the exact preprocessing pipeline, e.g., if using GLM residuals).
    >
    > **Practical Use:** Generate reliability-weighted maps:
    > \[ M^{\text{reliab}}_{q,v} = \rho_{q,v} \; \beta_q \, \Delta_{q,v} \]

    > **Context: Implementation Cost & Strategies for ρ_q,v**
    >
    > **Challenge:** Storing Δ_q,v for every fold (S), contrast (Q), voxel (V) in every searchlight (#SL) can consume significant RAM (e.g., S × Q × V × #SL × 8 bytes).
    >
    > **Cost-Saving Strategies:**
    > 1.  **On-the-fly Welford variance:** Update the running mean and sum of squared differences (M2) for Δ_q,v as each fold's Δ_fold arrives. Discard Δ_fold afterwards. Requires only O(Q × V) memory buffer per searchlight.
    >     ```R
    >     # Pseudo-R for Welford inside sphere loop
    >     # Needs mean_D, m2_D (both QxV) initialized
    >     delta <- Delta_fold - mean_D
    >     mean_D <- mean_D + delta / fold_id
    >     m2_D <- m2_D + delta * (Delta_fold - mean_D)
    >     # After loop: var_D <- m2_D / (S - 1)
    >     ```
    > 2.  **Center-voxel only:** Compute and store ρ_q,v only for the center voxel of the searchlight. Reliability information gets implicitly smoothed across the brain by overlapping searchlights.
    > 3.  **Split-half reliability (S=2):** If using only two folds (e.g., odd vs. even runs), variance is analytic: Var = (Δ₁ - Δ₂)² / 4. Matches Walther et al. (2016) approach. Coarser but very cheap.
    > 4.  **Jack-knife (Leave-One-Run-Out, S >= 3):** Recompute Δ leaving out run s (Δ_(-s)). Variance is computed from these S pseudo-values. Requires minimal extra buffering if Δ_(-s) is recomputed on demand.
    >     `jack_var = (S-1)/S * sum( (Delta_neg_s - mean(Delta_neg_s))^2 )`
    > 5.  **Bootstrap (Fold-level Resampling, S >= 3):** Sample S folds with replacement B times (e.g., B=100). Compute variance of Δ across bootstrap samples. CPU scales with B, but uses pre-computed per-fold Δ (if stored/cached) or recomputes on the fly.
    >
    > **Recommendation:** Welford or Jack-knife offer good precision with low memory overhead for S >= 3. Split-half is viable for S=2 or when coarse reliability is sufficient.

    <!-- FUTURE IMPLEMENTATION NOTES FOR RELIABILITY (ρ_q,v) -->
    <details>
    <summary><strong>[Future Plan] Detailed Implementation Notes for Reliability (ρ_q,v)</strong></summary>

    Based on the discussion on implementation efficiency and the blueprint provided:

    **Guiding Principles:**
    1.  Avoid modifying `compute_crossvalidated_means_sl` to return fold-wise data.
    2.  Compute fold-specific Δ (`Delta_fold_sl`) on-the-fly within `train_model.contrast_rsa_model`.
    3.  Use Welford's online algorithm for variance calculation to minimize memory.
    4.  Gate the feature behind a `calc_reliability` flag in the `contrast_rsa_model` constructor (defaulting to `FALSE`).

    **Proposed `rho` Formula:**
    Use the ICC-like form `rho = sigma2_noise_param / (var_delta + sigma2_noise_param)` where `sigma2_noise_param = (S-1) * var_delta` for `S > 1` (number of folds).
    This simplifies to `rho = (S-1)/S` for `var_delta > 0`. Handle edge cases:
    - If `var_delta = 0` (perfect stability), `rho = 1`.
    - If `S = 1`, `rho = 0` (unless `var_delta = 0`, then `rho = 1`).

    **Implementation Steps:**

    **1. `contrast_rsa_model` Constructor (`R/contrast_rsa_model.R`):**
        *   Add `calc_reliability = FALSE` parameter.
        *   Store in model spec via `create_model_spec`.
        *   Update Roxygen docs.

    **2. `train_model.contrast_rsa_model` (`R/contrast_rsa_model.R`):**
        *   Initialize `rho_vc_sl = rep(1, Q)` (default: neutral weight).
        *   **If `obj$calc_reliability` is `TRUE`:**
            *   Get `S = get_nfolds(cv_spec)`.
            *   Get `C`, `V_sl`, `Q`, `mvpa_des`.
            *   Initialize Welford accumulators: `m_Delta_sl = matrix(0, V_sl, Q)`, `M2_Delta_sl = matrix(0, V_sl, Q)`.
            *   Loop `s_idx` from 1 to `S`:
                *   Get `train_indices_fold_s`.
                *   Handle empty/problematic folds (skip or ensure NA propagation).
                *   Compute `U_hat_sl_fold_s` (K x V_sl) using logic similar to `compute_cv_means` (subsetting `sl_data`, `aggregate`, re-indexing, handling missing conditions -> NAs).
                *   If `anyNA(U_hat_sl_fold_s)`, skip Welford update for this fold or ensure NAs propagate correctly.
                *   `Delta_fold_sl = t(U_hat_sl_fold_s) %*% C` (V_sl x Q).
                *   Perform Welford update on `m_Delta_sl` and `M2_Delta_sl` using `Delta_fold_sl`.
            *   After loop (let `S_eff` be number of valid folds processed):
                *   Initialize `rho_sl = matrix(0, V_sl, Q)`.
                *   If `S_eff > 1`:
                    *   `var_Delta_sl = M2_Delta_sl / (S_eff - 1)`.
                    *   `sigma2_noise_param_sl = (S_eff - 1) * var_Delta_sl`.
                    *   `denominator_rho = var_Delta_sl + sigma2_noise_param_sl`.
                    *   `rho_sl = sigma2_noise_param_sl / denominator_rho`.
                    *   Handle `denominator_rho < 1e-10` (set `rho_sl = 1`).
                    *   Handle `NA` values in `rho_sl` (set `rho_sl = 0`).
                *   Else if `S_eff == 1`:
                    *   Handle perfect stability case: `rho_sl[M2_Delta_sl == 0] = 1`.
            *   Extract `rho_vc_sl = rho_sl[center_idx, , drop = TRUE]`.
        *   Modify final metric calculation:
            ```R
            final_delta_vc_sl = delta_vc_sl # Potentially normalized delta
            if (obj$output_metric == "beta_delta") {
                # Apply reliability weight
                result_vector <- beta_sl * final_delta_vc_sl * rho_vc_sl 
            } else if (obj$output_metric == "delta_only") {
                # Keep raw (potentially normalized) delta? Or apply rho?
                # Decision needed: For now, keep as is.
                result_vector <- final_delta_vc_sl 
            } else { # beta_only
                result_vector <- beta_sl
            }
            ```
    </details>
    <!-- END FUTURE IMPLEMENTATION NOTES -->

- [x] **Voxel-Specific RDM Reconstruction (r_v) (D.2):** Implement the calculation of the voxel-specific RDM reconstruction score (`r_v`) as the `"recon_score"` output metric.

    > **Context: Why Voxel Reconstruction Score (r_v)?**
    >
    > This metric answers: *How much does this single voxel v, with its specific profile of contributions across all contrasts Q, matter for reproducing the overall representational geometry (Ĝ_empirical) observed in this searchlight?*
    >
    > **Intuition:** It's a "leave-one-voxel-in" score. If r_v ≈ 1, this single voxel carries much of the information about pairwise condition distances relevant to the contrasts. If r_v ≈ 0, the geometry relies on the pattern across *other* voxels.
    >
    > **Building the One-Voxel Model RDM (Ĝ^(v)):**
    > Use the voxel's signed contributions (Δ_q,v) scaled by the overall contrast importance (β_q) for the searchlight. Project this Q-dimensional profile back into the KxK condition space using the contrast matrix C:
    > \[ \hat{G}^{(v)} = C \; \text{diag}(\beta_1 \Delta_{1,v}, \dots, \beta_Q \Delta_{Q,v}) \; C^T \]
    >
    > **Calculating the Score (r_v):**
    > Correlate the vectorized lower triangle of the one-voxel RDM (Ĝ^(v)) with the vectorized lower triangle of the empirical RDM (Ĝ_empirical) from the searchlight:
    > \[ r_v = \text{corr}\left( \text{vec}_{\text{lower}}(\hat{G}^{(v)}), \; \text{vec}_{\text{lower}}(\hat{G}_{\text{empirical}}) \right) \]
    >
    > **Interpretation:** A map of r_v (or |r_v|) highlights voxels that are individually highly informative about the multi-contrast representational structure.
    >
    > **Implementation Cost:** Very low. Requires C, β_q, Δ_q,v, and vec(Ĝ_empirical), most of which are already computed. The matrix multiplications per voxel are small (KxQ, QxQ, QxK). The correlation is between two vectors of length K(K-1)/2.
    > ```R
    > # Pseudo-R for r_v inside sphere loop, after beta & Delta computed
    > y_emp <- lower(G_hat_empirical) # Precompute
    > beta_diag <- diag(beta)         # QxQ
    > r_v <- numeric(V)
    > for (v in 1:V) {
    >    # Delta is V x Q or Q x V depending on implementation
    >    # Assuming Delta is V x Q here:
    >    voxel_loading_diag <- diag(beta * Delta[v, ]) # Corrected: element-wise product
    >    G_hat_v <- C %*% voxel_loading_diag %*% t(C)
    >    r_v[v] <- cor(lower(G_hat_v), y_emp)
    > }
    > ```

## Additional Output Maps (Section E of Proposal)

- [x] **Direction-Only Maps (~M_q,v) (E.1):** Add an `output_metric` option (`"beta_delta_norm"`) to `contrast_rsa_model` to generate maps based on normalized Δ_q (`β_q * ~Δ_{q,v}`).
- [ ] **Reliability-Weighted Maps (M_q,v_reliab) (E.1):** Add an `output_metric` option (e.g., `"beta_delta_reliable"`) to generate maps weighted by voxel reliability (`ρ_{q,v} * β_q * Δ_{q,v}`). Requires completion of D.1.
- [x] **Composite Map (w_v) (E.2):** Add an `output_metric` option (`"composite"`) to calculate a composite map representing the net pull of a voxel (e.g., `Σ_q (β_q * ~Δ_{q,v})`), including orthonormality check.

## Documentation & Examples

- [ ] **Add Examples (`@examples`):** Add runnable examples to all exported functions, particularly `contrast_rsa_model()` and `run_searchlight.contrast_rsa_model()`, creating necessary dummy objects.
- [ ] **Update Documentation:** Ensure all new parameters and options added during Phase 2 are clearly documented.

## Refactoring & Code Quality

- [x] **Refactor `wrap_out`/`create_searchlight_performance`:** Simplify the interaction between these functions in `R/searchlight.R`. (Done by making `wrap_out` directly create spatial objects and commenting out `create_searchlight_performance`, and adjusting `print.searchlight_result`).
- [x] **Check `cv_spec` Validity:** Implement the `TODO` check for `cv_spec` validity in `R/compute_cv_means.R`.
- [ ] **Verify `DESCRIPTION`:** Double-check that all necessary package dependencies (`neuroim2`, `neurosurf`, `dplyr`, `purrr`, etc.) are correctly listed in the `Imports` field of the `DESCRIPTION` file.
- [x] **Review Block Exclusion Logic (C4):** Addressed by deriving `condition_block_list` in `msreve_design` and using it in `train_model.contrast_rsa_model`.
- [x] **Replace `aggregate` with `rowsum` in `compute_cv_means` (Audit Quick Win #2):** Replaced `stats::aggregate` with `rowsum` and `table` for calculating fold-wise condition means for performance.

---

## Appendix: Context on Interaction Effects (Section F)

> **Why Interactions Matter in RSA:**
> Main effect contrasts test simple geometric predictions (e.g., animate vs. inanimate separation). Interactions test if the geometry has extra structure beyond the sum of main effects.
>
> **Construction:** An interaction contrast column (c_inter) is the element-wise product (⊙) of two centered main effect columns (c_p, c_q): `c_inter = c_p ⊙ c_q`.
>
> **Interaction RDM:** The model RDM for the interaction (`R_inter = c_inter %*% t(c_inter)`) predicts similarity (+1) for condition pairs sharing the same interaction sign (e.g., both are [+p, +q] or both are [-p, -q]) and dissimilarity (-1) for pairs differing in interaction sign (e.g., [+p, +q] vs. [+p, -q]). This captures a pattern orthogonal to both main effect RDMs.
>
> **Orthogonalization:** It's crucial to orthogonalize the full contrast matrix `C_exp = [C_main | C_inter]` (e.g., using `rMVPA::contrasts(..., orth=TRUE)` or `rMVPA::orthogonalize_contrasts()`). This ensures that the regression coefficient β_inter reflects variance uniquely explained by the interaction term after accounting for main effects.
>
> **Interpretation:**
> - `Δ_inter,v = t(U_hat) %*% c_tilde_inter`: Voxel v's projection onto the *orthogonalized* interaction contrast.
> - `β_inter * Δ_inter,v`: How strongly voxel v contributes *uniquely* to the interaction pattern.
> - Example: A voxel might respond strongly only to "animate-small" conditions, contributing positively to an animacy-by-size interaction, even if its response to animacy or size alone isn't distinct.
>
> **Workflow:** Define interactions in the `spec` for `rMVPA::contrasts()`, ensure `orth=TRUE`, then run `contrast_rsa_model`. The engine handles interaction columns like any other contrast. 