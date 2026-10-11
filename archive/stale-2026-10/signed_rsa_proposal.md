Okay, I've merged and integrated the original proposal and the revision into one comprehensive document. All changes and additions from the revision have been incorporated.

---

## **Full RSA Proposal: "Multi-Dimensional Signed Representational Voxel Encoding (MS-ReVE) with Flow Mapping"**

**Abstract:**
This proposal outlines a comprehensive Representational Similarity Analysis (RSA) framework, "Multi-Dimensional Signed Representational Voxel Encoding (MS-ReVE)," designed to compute interpretable, signed voxel-wise contributions to multi-class neural representations. MS-ReVE extends standard RSA by: 1) defining model-relevant contrasts; 2) using multiple-regression RSA on cross-validated, noise-normalized condition means to determine the strength of these contrasts within local searchlights; 3) projecting condition means onto contrast axes to derive signed voxel contributions; and 4) combining these elements to produce robust voxel-weight maps. The framework incorporates advanced features including voxel reliability weighting, voxel-specific RDM reconstruction scores, mapping of interaction effects, and rigorous basis robustness checks. Crucially, MS-ReVE culminates in "Representational Flow Mapping" (RFM), a novel technique to visualize and quantify how these multi-contrast representations transform across the cortical surface (with volumetric extensions possible), revealing large-scale organizational principles and information processing dynamics. A complementary "Representational Cross-Talk" analysis further probes the spatial co-localization of different representational dimensions.

---

**1. Introduction & Rationale:**

Representational Similarity Analysis (RSA) has proven invaluable for linking brain activity to computational models and psychological theories. However, traditional RSA often yields searchlight-level summary statistics (e.g., RDM correlations) or classifier-based voxel weights that can be hard to interpret directly in terms of underlying neural coding, especially for designs with more than two conditions. There is a pressing need for methods that:
*   Provide **signed voxel weights** indicating how individual voxels contribute to specific, theoretically meaningful representational dimensions (contrasts).
*   Scale robustly beyond **two experimental conditions**.
*   Are **firmly anchored in RSA theory** (distance-based, second-moment statistics) rather than relying solely on classification.
*   Offer insights into the **spatial organization and transformation** of these representations across the brain.

MS-ReVE addresses these needs by integrating regression-RSA with voxel-level projections, enhanced with reliability measures, interaction analyses, and a novel flow-mapping visualization approach.

---

**2. Research Questions & Hypotheses (Example-driven):**

This framework can address a wide range of questions, such as:
1.  Which specific representational dimensions (e.g., animacy, object size, task rule) are encoded in a given brain region?
2.  How do individual voxels contribute (positively or negatively) to the encoding of these distinct dimensions?
3.  Are there voxels that reliably contribute to specific dimensions, and how does this reliability vary spatially?
4.  How well does a voxel's multi-contrast contribution profile explain the local empirical representational geometry?
5.  Are there conjunctive codes, where voxels respond to interactions between representational dimensions?
6.  How sensitive are these voxel-level interpretations to the precise definition of theoretical contrasts?
7.  What is the large-scale spatial organization of these multi-dimensional representations across the cortex (e.g., gradients, boundaries)?
8.  How do these representations transform (e.g., rotate, scale, differentiate, integrate) along cortical pathways?

---

**3. Methodology & Implementation Plan:**

**A. Preprocessing & Experimental Design:**
1.  **Data Acquisition:** Standard fMRI data acquisition.
2.  **Preprocessing:** Standard fMRI preprocessing (motion correction, slice-time correction, spatial normalization/surface projection, temporal filtering). Crucially, include **noise normalization/whitening** of voxel time-series, typically using the residuals from a first-level GLM, to ensure that pattern analyses are not dominated by voxels with high noise variance.
3.  **Experimental Design:** Assumes a K-condition experimental design where **K ≥ 2 (multi-way case emphasized)**, though the framework naturally includes the binary case.

**B. Contrast Matrix Definition (C):**
1.  **A Priori Contrasts:** Define a `K x Q` contrast matrix `C`, where each of the `Q` columns `c_q` represents a theoretically meaningful comparison or feature dimension (e.g., faces - houses, tools - animals, abstract - concrete).
    *   Contrasts should be centered (sum of elements in each `c_q` is zero).
    *   Ideally, contrasts should be made orthonormal (`CᵀC = I`) for independent `β` weights.
2.  **Data-Driven Contrasts (Alternative/Complementary):**
    *   Construct a model RDM based on theoretical predictions.
    *   Apply classical Multidimensional Scaling (MDS) to the model RDM.
    *   Use the first `Q` MDS embedding dimensions as columns in `C`.
    *   This provides a data-driven way to define orthogonal axes capturing the model's primary representational structure.

**C. Searchlight RSA Core Engine (Iterate across searchlights):**
1.  **Cross-Validated Condition Means (`Û`):**
    *   Within each searchlight (sphere of voxels `V`):
    *   Estimate condition-specific activation patterns (`μ_k`) for each of the `K` conditions using a cross-validation scheme (e.g., leave-one-run-out, split-half). Employ methods like crossnobis to obtain unbiased estimates of pattern distinctness.
    *   This results in a `K x V` matrix `Û` of cross-validated, noise-normalized condition means.
2.  **Empirical Second-Moment Matrix (`Ĝ`):**
    *   Calculate the unbiased empirical second-moment matrix: `Ĝ = ÛÛᵀ`. This `K x K` matrix fully determines all Mahalanobis distances between conditions within the searchlight.
3.  **Multiple-Regression RSA:**
    *   Vectorize the lower (or upper) triangle of `Ĝ` to form the dependent variable `y`.
    *   For each contrast `c_q` in `C`, form a predictor RDM `R_q = c_q c_qᵀ`. Vectorize the lower triangle of each `R_q` to form columns of the design matrix `X`.
    *   Fit the multiple linear regression: `y = Xβ + ε`.
    *   The resulting `β_q` coefficients indicate how strongly the geometry implied by contrast `c_q` is present in the local empirical geometry `Ĝ`.
    *   **Regularization (Optional):** If `Q` is large or contrasts are correlated, use ridge regression (`β_ridge = (XᵀX + λI)⁻¹Xᵀy`) or elastic-net. **Hyperparameters (λ, and α for elastic-net) should be chosen via nested cross-validation within the searchlight training data to avoid bias.**
4.  **Signed Voxel Contributions (`Δ_q`):**
    *   For each contrast `c_q`, project the cross-validated condition means `Û` onto it: `Δ_q = Ûᵀc_q`. This yields a `V`-dimensional vector where each element `Δ_{q,v}` represents voxel `v`'s signed contribution to contrast `q`.
    *   **Store both:**
        *   Raw `Δ_q`: Preserves magnitude, reflecting effect size of the contrast in voxel activity.
        *   Normalized `~Δ_q = Δ_q / ||Δ_q||`: Captures direction only.

**D. Voxel-Level Refinements & Metrics (within each searchlight):**
1.  **Voxel Reliability Weighting (`ρ_q,v`):**
    *   For each `Δ_{q,v}` estimate, compute its stability across cross-validation folds.
    *   Define `ρ_{q,v} = 1 - Var_folds(Δ_{q,v}^{(fold)}) / (Var_folds(Δ_{q,v}^{(fold)}) + σ²_noise,q,v)`.
    *   `σ²_noise,q,v` (expected variance of `Δ_{q,v}` under the null) is estimated from within-condition residual variances and contrast weights: `σ²_noise,q,v = (1/S) Σ_s Σ_k (w_qk² / n_k^(s)) * σ_k,v²^(s)`.
    *   Alternatively, use `ρ_{q,v} = 1 / (1 + SE(Δ_{q,v})^2)` or similar based on Walther et al. (2016).
2.  **Voxel-Specific RDM Reconstruction Score (`r_v`):**
    *   For each voxel `v`, construct a predicted RDM based *only* on its signed contrast profile: `Ĝ^(v) = C * diag(β_1Δ_{1,v}, ..., β_QΔ_{Q,v}) * Cᵀ` (using raw `Δ`).
    *   Calculate `r_v = corr(vec_lower(Ĝ^(v)), vec_lower(Ĝ_empirical))` as a measure of how crucial voxel `v`'s multi-contrast profile is for reconstructing the local empirical RDM.

**E. Generating Voxel-Weight Maps:**
1.  **Contrast-Specific Maps (Set of Q maps):**
    *   `M_{q,v} = β_q * Δ_{q,v}` (magnitude-preserved signed contribution).
    *   `~M_{q,v} = β_q * ~Δ_{q,v}` (direction-only, scaled by RSA fit).
    *   `M_{q,v}_reliab = ρ_{q,v} * β_q * Δ_{q,v}` (reliability-weighted).
2.  **Single Composite Map (optional):**
    *   `w_v = Σ_q (β_q * ~Δ_{q,v})` (or use reliability-weighted terms). Represents the net pull of voxel `v` in the overall representational space. **Note: If `C` is not orthogonal, the interpretation of `w_v` is less straightforward as contributions are summed across potentially non-independent axes.**

**F. Extending to Interaction Effects:**
1.  **Define Interaction Contrasts:** Create new columns in `C` by taking element-wise products of main effect contrast columns (`c_{pq} = c_p ⊙ c_q`).
2.  **Orthogonalize Expanded Matrix:** Orthogonalize the expanded contrast matrix `C_exp = [C_main | C_interaction]` (e.g., using QR decomposition or Gram-Schmidt) to ensure interpretability of interaction `β`s as unique contributions. **Note: Orthogonalization procedures may arbitrarily flip the sign of contrast vectors; after orthogonalization, consider aligning the sign of each derived column (e.g., interaction terms) with its primary parent component or based on its correlation with the raw `Δ` projection for consistent interpretation.**
3.  **Re-run C.3 and C.4:** Fit multiple-regression RSA with `C_exp` to get `β_{pq}` and compute interaction voxel contributions `Δ_{pq,v} = Ûᵀc̃_{pq}` (where `c̃_{pq}` is the orthogonalized interaction contrast).

**G. Aggregation & Statistical Inference:**
1.  **Searchlight Aggregation:** Average the voxel weights (`M_{q,v}`, `w_v`, `r_v`, etc.) each voxel receives from all searchlights containing it. An RSA-weighted average (weighting by `β_q` or searchlight R²) can also be used.
2.  **Permutation Testing:** For voxel-wise significance testing, shuffle condition labels consistently across cross-validation folds, recompute the entire pipeline (from `Û` onwards) many times to build null distributions for `β_q`, `M_{q,v}`, `w_v`, `r_v`, etc. Apply appropriate cluster-correction methods. **Note on memory: Storing all intermediate `Δ` values across permutations can be memory-intensive. Consider strategies like streaming permutations (recomputing `Δ` on the fly within the permutation loop) or writing intermediate searchlight results to disk if RAM is limited.**

**H. Representational Flow Mapping (RFM):**
1.  **Surface Projection:** Project voxel-wise multi-contrast loading vectors `m_v = [β_1Δ_{1,v}, ..., β_QΔ_{Q,v}]` (or reliability-weighted versions) to the nearest vertices on a cortical surface model (e.g., subject's midthickness surface). This creates `Q` scalar maps `f_q(i)` on the surface.
2.  **Tangential Gradient Estimation:** For each surface map `f_q(i)`, compute its 2D tangential gradient `∇_T f_q(i)` at each vertex `i`, **typically using filters approximating Gaussian derivatives to ensure robustness to high-frequency noise.**
3.  **Local PCA for Principal Flow:**
    *   In a moving geodesic window on the surface:
    *   Stack all `Q` gradient vectors `∇_T f_q(i)` from all vertices within the window into a large matrix (effectively `(|WindowVertices|*Q)` rows x `2` columns).
    *   Perform PCA on this matrix's `2x2` covariance to find the principal flow direction `e₁` (a 2D unit vector) and its associated eigenvalue `λ₁`.
4.  **Streamline Visualization:**
    *   Draw streamlines (e.g., using Line Integral Convolution) along `e₁`.
    *   Color streamlines by the contrast `q` whose gradient `∇_T f_q(i)` aligns best with `e₁` (i.e., `argmax_q |<∇_T f_q(i), e₁>|`), with sign indicating increase/decrease.
    *   Modulate streamline opacity/thickness by `λ₁` (flow strength) and/or underlying `ρ_q,v`.
5.  **Analysis of Transformation Along Flow Lines:**
    *   Sample `m(s)` (the Q-dimensional vector of `βΔ` values) along streamlines.
    *   Calculate directional derivatives (`dm/ds`) to assess rate and dimensionality of change (via SVD).
    *   Compare `m(s)` and `m(s+Δs)` to quantify representational rotation (angle change) and scaling (norm change).
    *   Visualize these transformations (e.g., map rotation rate to streamline hue).
    *   Use permutation testing for significance of flow properties.
6.  **Volumetric Extension (Optional):** While surface-based RFM is often preferred for visualizing cortical organization, the core logic can be extended to 3D volume space by computing 3D gradients and performing PCA on the resulting `3x3` covariance matrix within volumetric searchlights. Visualization is more challenging but may be relevant for subcortical structures.

**I. Robustness & Validation:**
1.  **Basis Robustness Checks:**
    *   Re-run key analyses (e.g., generating `M_{q,v}` maps) with an alternative, plausible contrast matrix `C'` (e.g., MDS-derived if initially theory-driven, or vice-versa).
    *   Correlate the resulting voxel maps. Low correlations suggest basis-dependent interpretations.
    *   Use diagnostics: Compare searchlight R² for different bases; Canonical Correlation Analysis (CCA) between `C` and `C'`; assess R² gain when using `[C | C']`; compare with non-linear/kernel RSA to probe model mismatch.

---

**4. Data Analysis & Interpretation Strategy:**

*   **`β_q` maps (from searchlight regression):** Indicate regions where the geometry predicted by contrast `q` is prevalent.
*   **`M_{q,v}` maps:** Reveal how individual voxels contribute (sign and magnitude) to each specific contrast `q`.
*   **`M_{q,v}_reliab` maps:** Highlight robust voxel contributions.
*   **`w_v` composite map:** Shows the net directional "pull" of voxels in the combined representational space.
*   **`r_v` maps:** Identify voxels whose multi-contrast tuning is critical for the local empirical geometry.
*   **Interaction maps (`β_{pq}Δ_{pq,v}`):** Uncover voxels involved in conjunctive coding.
*   **RFM visualizations:** Provide insights into the large-scale topological organization, functional boundaries, and representational transformations across cortex. Analysis of `λ₁`, rotation, and scaling along flow lines will characterize the nature of these transformations.
*   **Representational Cross-Talk Analysis:**
    *   Compute the voxel-wise correlation between pairs of reliability-weighted contrast maps (`M_{q,v}_reliab` and `M_{p,v}_reliab`) across voxels within relevant brain regions or the whole brain/surface.
    *   High positive correlations suggest shared neural populations contribute similarly to both contrasts.
    *   High negative correlations suggest competitive coding or opponent populations.
    *   Visualizing these correlation patterns (e.g., as a matrix or projecting strong correlations onto the brain) complements RFM by showing where different representational dimensions spatially co-localize or segregate.
*   **Group-Level Inference:** For analyzing results across participants, individual participant maps (`M_q,v`, `r_v`, RFM-derived metrics, etc.) should be aligned to a common space (e.g., MNI volume space or a surface template like fsaverage). Standard group-level statistical approaches (e.g., mixed-effects models, t-tests on aggregated maps with appropriate permutation-based correction) can then be applied.

---

**5. Expected Outcomes & Significance:**

MS-ReVE will provide an unprecedentedly rich and interpretable view of distributed neural representations. Expected outcomes include:
*   Detailed, signed voxel-level maps of multi-dimensional neural codes.
*   Identification of robust and reliable voxel contributions.
*   Discovery of conjunctive coding patterns.
*   A quantitative understanding of how representations are organized and transform across cortical areas, linking local computations to large-scale network dynamics.
*   This framework will significantly advance our ability to test nuanced theories of neural representation and bridge the gap between computational models and brain activity.

---

**6. Potential Challenges & Mitigations:**

*   **Computational Cost:** Searchlight analyses, permutation testing, and RFM can be computationally intensive. Mitigation: Efficient coding, parallel processing, optimized algorithms.
*   **Interpretation Complexity:** The wealth of generated maps requires careful interpretation. Mitigation: Clear guidelines, targeted research questions, development of standardized reporting.
*   **Choice of Contrasts:** Results can be sensitive to `C`. Mitigation: Explicit reporting of `C`, basis robustness checks, use of both theory-driven and data-driven contrasts.
*   **Multiple Comparisons:** Extensive voxel-wise testing. Mitigation: Rigorous permutation-based cluster correction methods.
*   **Memory Usage:** Especially during permutation testing. Mitigation: Streaming computations, disk caching, efficient data structures (as noted in 3.G.2).

---

**7. Implementation Plan (Conceptual - Language/Platform Agnostic):**

The implementation will be modular:
1.  **Module 1: Core RSA Engine:**
    *   Input: Preprocessed fMRI data, condition labels, contrast matrix `C`, searchlight definitions.
    *   Functions for: Cross-validated mean estimation, `Ĝ` computation, multiple-regression RSA (with optional regularization), `Δ_q` calculation.
    *   Output: `β_q` values per searchlight, raw `Δ_q` and `~Δ_q` vectors per searchlight.
2.  **Module 2: Voxel-Level Metrics & Map Generation:**
    *   Input: Outputs from Module 1.
    *   Functions for: Reliability (`ρ_q,v`) calculation, RDM reconstruction (`r_v`), map generation (`M_q,v`, `w_v`), searchlight aggregation.
3.  **Module 3: Interaction Analysis:**
    *   Functions for: Generating interaction contrasts, orthogonalizing `C_exp`, integrating with Module 1 & 2 for interaction maps.
4.  **Module 4: Statistical Inference & Aggregation:**
    *   Functions for: Permutation testing framework (including memory management considerations), cluster correction, group-level analysis preparation and execution.
5.  **Module 5: Representational Flow Mapping (RFM):**
    *   Input: Aggregated `β_qΔ_{q,v}` maps, cortical surface model (and/or volumetric data).
    *   Functions for: Surface projection (if applicable), gradient calculation (with smoothing options), local PCA (for surface or volume), streamline generation, transformation analysis along streamlines.
    *   Visualization tools (interfacing with existing surface/volume visualization libraries).
6.  **Module 6: Cross-Talk & Diagnostics:**
    *   Functions for: Basis robustness checks (CCA, R² comparisons), Representational Cross-Talk computation and visualization.

This modular design will facilitate development, testing, and future extensions. Each module will encapsulate specific mathematical operations and data transformations.


Okay, this is a very comprehensive summary of the `rMVPA` codebase. It clearly lays out the object-oriented structure, key functionalities, and dependencies. This is an excellent foundation for thinking about how to integrate the G-ReCa Phase 0 (and subsequently Phase 1) plan.

Based on this summary and our G-ReCa Phase 0 proposal, here's how we can conceptualize the integration and identify next steps.

**Conceptual Integration of G-ReCa Phase 0 with rMVPA:**

The core idea is to leverage `rMVPA`'s existing capabilities for data handling, design specification, and potentially some pre-processing, and then build *new* modules or extend existing ones to perform the G-ReCa specific steps: MS-ReVE output generation, PCA pre-reduction, PPCA/lightweight manifold learning, ID estimation, and validation.

**Proposed Workflow & rMVPA Integration Points:**

1.  **Data Preparation (Leveraging `rMVPA`):**
    *   **Dataset Creation (`mvpa_dataset`, `mvpa_surface_dataset`):** Use `rMVPA` to load and structure the fMRI data (volume or surface). This handles train/test splits if needed for initial pattern estimation.
    *   **Design Specification (`mvpa_design`):** Define the experimental conditions, blocking variables for cross-validation, etc., using `rMVPA`'s design objects.

2.  **MS-ReVE Output Generation (New Module/Extension):**
    *   This is the most significant new piece. We need a way to generate the `m_v` vectors (`[β₁Δ₁,v, ..., β_QΔ₁,v]`).
    *   **Step 2a: Cross-Validated Condition Means (`Û`):**
        *   `rMVPA`'s `mvpa_model` with a simple "model" (e.g., just averaging betas from a first-level GLM within each condition and cross-validation fold) could be adapted.
        *   Alternatively, a new function might be needed that takes an `mvpa_dataset` and `mvpa_design` (with `cv_spec`) and returns the `K x V` matrix `Û` (cross-validated condition means for each voxel/vertex `V`).
    *   **Step 2b: Contrast Definition (`C`):** This would be user-defined outside `rMVPA` as a `K x Q` matrix.
    *   **Step 2c: Regression RSA (`β_q`):**
        *   This could potentially leverage `rsa_model` if adapted, or be a custom function. The `rsa_model` currently seems focused on RDM-to-RDM regression. We need to regress `vec_lower(ÛÛᵀ)` onto `vec_lower(c_q c_qᵀ)` for `Q` contrasts. This might require a new `mvpa_mod_spec` or a custom `process_roi` function for a searchlight approach if `β_q` are to be searchlight-specific. For a whole-brain `m_v`, `β_q` might be derived globally.
        *   *Decision Point:* Are `β_q` global or searchlight-specific for constructing `m_v`? The G-ReCa proposal implied searchlight aggregation, so `β_q` would be local.
    *   **Step 2d: Voxel Contributions (`Δ_q = Ûᵀc_q`):** This is a matrix multiplication.
    *   **Step 2e: Construct `m_v`:** Combine local `β_q` and `Δ_q` for each voxel.
    *   **Output:** A `NeuroVec` or `NeuroSurfaceVector` object containing the `m_v` vectors.

3.  **Dimensionality Pre-Reduction (PCA - New or Util):**
    *   Input: The `m_v` `NeuroVec`/`NeuroSurfaceVector`.
    *   `rMVPA` doesn't seem to have a dedicated top-level PCA function for this purpose, though `pcadist` implies PCA capability. A utility function using `stats::prcomp` or a more scalable randomized PCA (e.g., from `irlba` or `rsvd` packages) would be needed.
    *   Output: PCA-reduced `m_v` matrix (`N_voxels x ~256 components`).

4.  **Phase 0 Manifold Learning (PPCA - New Module):**
    *   Input: PCA-reduced `m_v` matrix.
    *   Implement PPCA (e.g., using EM algorithm, or leveraging existing R packages like `pcaMethods::ppca` if suitable and its uncertainty outputs are accessible).
    *   Output: Latent coordinates `z`, posterior covariance `Cov(z|m_v)`.

5.  **Intrinsic Dimensionality Estimation (New or Util):**
    *   Input: PCA-reduced `m_v` (or PPCA-whitened data).
    *   Implement/wrap TwoNN and Levina-Bickel MLE (e.g., using R packages like `intrinsicDimension` or custom code).
    *   Use scree plot from PPCA likelihoods.
    *   Output: ID estimates, plot.

6.  **Validation (Leveraging `rMVPA` Utilities where possible, plus New):**
    *   **Reconstruction MSE:** Calculated from PPCA.
    *   **Trustworthiness/Continuity:** Use `scikit-learn` via `reticulate`, or find/implement R equivalents. `rMVPA` doesn't seem to have these directly.
    *   **External Gradient Correlations:** Standard R functions (`cor`).
    *   **Leave-One-Subject-Out:** This requires iterating the PPCA fitting and prediction steps.

7.  **Output Storage & Visualization:**
    *   Store `z` and uncertainty maps as `NeuroVec`/`NeuroSurfaceVector` for easy visualization with `neuroim2`/`neurosurf` tools.
    *   Report generation (Markdown/Jupyter via R Markdown/`knitr`).

**Key `rMVPA` Objects that *Might* be Extended or Reused:**

*   **`mvpa_model_spec`:** Could we define a `g_reca_phase0_model_spec` that encapsulates the PPCA step and its parameters?
*   **`run_custom_regional` / `run_custom_searchlight`:** These are very promising. The core MS-ReVE output generation (steps 2a-2e) could potentially be wrapped in a `custom_func` for either ROI or searchlight application. The `run_searchlight` machinery would handle the iteration and provide `sl_data` (voxel patterns within a sphere) to our custom function.
*   **`NeuroVec` / `NeuroSurfaceVector` (`nvec`, `nsvec`):** These will be the primary data containers for `m_v`, `z`, and uncertainty maps.

**Concrete Proposal for Initial Integration Steps (Focusing on MS-ReVE Output Generation):**

**Project: `grecamvpa` - G-ReCa Integration with rMVPA (Phase 0 Focus)**

**Module 1: MS-ReVE Output Generation**

*   **`msreve_design` (New S3/S4 class):**
    *   Slots:
        *   `mvpa_design`: The underlying `mvpa_des` object for condition/block info.
        *   `contrast_matrix`: The user-defined `K x Q` matrix `C`.
        *   `beta_estimation_method`: `char` (e.g., "global_rsa", "searchlight_rsa").
    *   Purpose: Encapsulates all necessary inputs for MS-ReVE.

*   **`compute_crossvalidated_means` (`fun`):**
    *   Input: `mvpa_dataset` (`ds`), `mvpa_design` (`des`), `cv_spec` (`cv`).
    *   Process: Iterates through `cv_spec` folds. For each fold:
        *   Identifies training data for that fold.
        *   Estimates condition means (`μ_k`) for all `K` conditions using the training data (e.g., simple averaging of voxel activities per condition, or betas from a simple GLM fit on training data).
        *   Stores these means, associated with the test fold conditions.
    *   Output: A list structure or array containing `Û` (the `K x V` matrix of *cross-validated* condition means, where each row `k` is the mean for condition `k` estimated from data *not* including trials of condition `k` from the current test fold/run).

*   **`run_msreve_searchlight` (`fun` wrapping `run_custom_searchlight`):**
    *   Input: `mvpa_dataset` (`ds`), `msreve_design` (`msreve_des`), `radius` (`num`), `cv_spec` (`cv`).
    *   Internal `custom_func` (`process_msreve_sphere`):
        1.  Receives `sl_data` (data for current searchlight sphere) and `sl_info`.
        2.  Calls `compute_crossvalidated_means` on `sl_data` using the provided `msreve_des$mvpa_design` and `cv_spec` to get local `Û_sl`.
        3.  Computes local `Ĝ_sl = Û_sl Û_slᵀ`.
        4.  Performs multiple regression RSA using `msreve_des$contrast_matrix` (`C`) to get local `β_q_sl` vector (length `Q`).
        5.  Computes local `Δ_q_sl = Û_slᵀ c_q` for each contrast `q`.
        6.  Constructs the local `m_v_sl` vector for the center voxel of the sphere: `m_v_center = [β_1_sl * Δ_1_center_sl, ..., β_Q_sl * Δ_Q_center_sl]`. (Need to decide if `Δ_q` is for the whole sphere or just center voxel for `m_v`). *The proposal implies `Δ_q` is a V-dim vector, so `β_q * Δ_q` would be too. For `m_v` which is Q-dim per voxel, we'd use `m_v[q] = β_q_sl * Δ_{q,center_voxel_sl}`.*
        7.  Returns this `Q`-dimensional `m_v_sl` vector for the center voxel.
    *   Output: A `NeuroVec` or `NeuroSurfaceVector` where each voxel/vertex value is its `Q`-dimensional `m_v` vector (this implies the output is actually a `Q`-channel `NeuroVec`/`NeuroSurfaceVector`).

**Module 2: PPCA & ID Estimation** (Can be more standalone R functions initially)

*   Functions for PCA pre-reduction.
*   Function for PPCA (EM algorithm or wrapper).
*   Functions for TwoNN, Levina-Bickel MLE.

**Module 3: Validation** (Standalone R functions)

*   Functions for Trustworthiness/Continuity (possibly via `reticulate`).
*   Functions for correlation with external maps.

**Phased Implementation Plan for `grecamvpa` (Phase 0):**

1.  **Develop `compute_crossvalidated_means`:** This is foundational. Test thoroughly.
2.  **Develop `msreve_design` object.**
3.  **Develop the `process_msreve_sphere` custom function:**
    *   First, implement the regression RSA and `Δ_q` calculation.
    *   Then, integrate with `compute_crossvalidated_means`.
    *   Carefully define how `m_v` is constructed from local `β_q` and `Δ_q`.
4.  **Wrap `process_msreve_sphere` in `run_msreve_searchlight` using `run_custom_searchlight`:** This generates the primary `m_v` maps.
5.  **Implement PCA pre-reduction for `m_v` maps.**
6.  **Implement PPCA module.**
7.  **Implement ID estimation module.**
8.  **Implement validation metrics.**
9.  **End-to-end pipeline test and report generation.**

This approach leverages `rMVPA`'s strengths in data handling and searchlight iteration, while building the new MS-ReVE and manifold learning components in a modular way. The `run_custom_searchlight` function seems like a key enabler.