# Pattern-first spatial reduced-rank MVPA (`pattern_model`) — assessment and implementation plan

**Status:** Phases 0-6 implemented within the boundaries below; local package and integration gates passed as recorded below; dependent PR review and hosted gates pending. Supported-rank tests (5b), support envelope/local noise (3b), and joint hierarchical fitting remain extensions.
**Date:** 2026-09-06
**Base commit:** `acccd31` (master)
**Scope:** assess the "pattern-first spatial reduced-rank MVPA" proposal and turn it into a staged, verifiable implementation plan for rMVPA.

---

## Part A. Assessment

### A.1 Verdict

The method is coherent and worth building. The two identities the proposal leans on were re-derived and hold:

- **Haufe identity.** With `W = Psi^{-1} A G^{-1}`, `z = W'x = t + W'eps`, `Cov(z) = Phi + G^{-1}`, `Cov(x,z) = A(Phi + G^{-1})`, so `Cov(x,z) Cov(z)^{-1} = A`. Requires `Cov(t, eps) = 0` and the covariance model, not Gaussianity.
- **Leave-one-voxel-out information.** `(Psi_{-v,-v})^{-1}` is the Schur complement of `P = Psi^{-1}`, which gives `G_{-v} = G - h_v h_v'` with `h_v = (PA)_{v:}' / sqrt(P_vv)`. The matrix determinant lemma then gives `I(t; x_v | x_{-v}) = -1/2 log(1 - h_v' K h_v)`, `K = (Phi^{-1} + G)^{-1}`. Note this is information about the task scores `t = C'y`, not `y` itself when `q > r`. `P_vv` is cheap under `Psi = D + UU'` via Woodbury.

The integration architecture in the proposal (a `pattern_model` analysis family on the `fit_roi()` / `roi_result()` / `output_schema()` contract, a thin `run_global.pattern_model()`, a pure numerical core, a feature-aligned graph object, an explicit target accessor, matrix-free Haufe) is the right shape for this package and was checked against the source. The proposal's specific claims about existing code are accurate (see A.3).

### A.2 Substantive corrections and risks

These change the plan; they do not invalidate the design.

1. **The noise model determines what "local prediction" can mean.** With `Psi = D + UU'` and `h` on the order of 10–20 global components, the decoder `Psi^{-1}A` can cancel only a handful of global noise modes. Restricting to an ROI then gives `Psi_RR ≈ D_R + (a few global modes)`, so the restricted regional predictor is close to naive-Bayes-within-ROI. A searchlight LDA with shrinkage exploits *local* noise correlations that this model cannot represent. This is the single most likely way the method loses to searchlights on decoding accuracy. Plan response: ship diag+low-rank first (it is what makes everything tractable), but (a) benchmark explicitly against shrinkage-LDA searchlights and thresholded PLS, (b) reserve a slot in the noise abstraction for a graph-local precision term (neighbour regression / sparse precision on the adjacency graph) as a Phase 3b option, and (c) state the limitation in the docs from day one.

2. **Tuning burden is the practical blocker.** Rank × sparsity × support-TV × signed-TV × noise rank inside nested blocked CV is a five-dimensional grid. v1 must collapse this: fix noise rank by a residual-spectrum heuristic (not CV), tie smoothing to sparsity by a fixed ratio with a small set of ratio options, and tune only a warm-started (rank, lambda) path on inner blocked folds. Rough cost per fit at n=500, p=100k, r=4 is ~2e8 flops per gradient step, so a 30-iteration warm-started solve is seconds in base R; a full outer×inner×path grid is minutes to tens of minutes. Measure before writing C.

3. **TV and the support envelope are Phase 3b, not v1.** Group-lasso row sparsity plus a graph-Laplacian quadratic smoothness penalty on `A` gives a smooth-plus-simple-prox objective that FISTA handles matrix-free with `Matrix` sparse products, in R, today. The `||A_v|| ≤ g_v` envelope with TV on `g` is a constrained non-smooth problem needing a proper primal-dual solver and almost certainly C. The scientific point (coherent support with signed fine-scale loadings) survives the substitution: the Laplacian penalty on the envelope can be added later without changing the API, since `penalty = list(sparse, support_smooth, signed_smooth)` names the *roles*, not the functional form. The existing `spacenet_tvl1` PDHG solver and edge builder (`R/classifiers.R:1050,1164`) are single-column TV-L1 and can be generalized to row-group TV when 3b arrives.

4. **Rotation pays off only when `r > 1`.** For a binary contrast `r = 1` and rotation is moot. Varimax views should be built after classification/regression are verified (Phase 4), orthogonal only in v1, oblique deferred.

5. **Generative vs discriminative fit.** The objective regresses `X` on `YC` (encoding direction). Under misspecification this need not minimize decoding loss. The proposal already says to select hyperparameters by held-out *decoding* loss; the plan makes that the only selection criterion and adds a discriminative-refit-of-the-head option (refit the small `r`-dimensional prediction head by logistic/ridge regression on `u`) as a cheap hedge that keeps `A` fixed.

6. **Multi-response metrics must be defined once.** Per-response predictive R² against the training-mean baseline, plus a chosen aggregate (mean over responses, optionally set-weighted for `feature_sets_design`), plus mean correlation as a secondary. Never report squared correlation as R².

7. **Sequential rank testing is a research problem, not a v1 feature.** Ship `rank_selected` (CV) and `component_stability` (split-half / fold agreement of subspaces via principal angles) first. `rank_supported` requires a validated higher-rank null and belongs in Phase 5b.

### A.3 What the codebase already provides (verified at `acccd31`)

| Claim in proposal | Verified | Evidence |
|---|---|---|
| `fit_roi()`/`roi_result()`/`output_schema()` contract is the modern plugin path | Yes | `R/allgeneric.R:231,257,381`; `R/data_roi_result.R:37`; 17 registered `fit_roi` methods |
| Only `run_global.mvpa_model` exists; no default | Yes | `R/global_analysis.R:236`; TODO at `R/allgeneric.R:215-219,1513-1516` |
| `run_global` checks `length(y) == N` (fails for matrix targets) | Yes | `R/global_analysis.R:260-261` |
| `cov(X)` dense p×p in global Haufe path | Yes | `R/global_analysis.R:321`; also `R/importance.R:128,142,156` |
| `mvpa_design` separates `cv_labels` from `targets`; `y_train()` returns `cv_labels` | Yes | `R/design.R:34,188-227,300`; `targets` is unvalidated and may be a matrix |
| No `targets_test` field | Yes | only `y_test` (`R/design.R:303`) |
| `feature_sets_design` carries `X_train`, `X_test`, block vars, `time_series`, row weights | Yes | `R/feature_sets_design.R:83-163`; `row_weights` in `R/feature_sets.R:151` |
| `tune_model` uses row-wise `rsample::bootstraps` | Yes | `R/model_fit.R:48-76` |
| duplicate-column screening and `spacenet_tvl1` special case in `train_model` | Yes | `R/model_fit.R:498-499,569` |
| `region_importance` is random-subset ablation | Yes | `R/importance.R:262` |
| `spatial_nmf` has graph machinery | Yes | `.prepare_graph` `R/spatial_nmf.R:327`; `build_voxel_adjacency` (6/18/26-nbr) `:368`; `build_graph_laplacian` `:436`; all unexported |
| `mvpa_config(mode = "global", model_spec = )` routes to `run_global` | Yes | `R/workflow_api.R:14,39-52,391` |
| `src/` has C registration | Yes | `src/init.c`, 3 routines in `dual_lda_cholupdate.c` |

**What the proposal missed (and the plan now uses):**

- **Existing reduced-rank precedent inside the package.** `remap_rrr_model` (`R/remap_rrr_model.R:80`, RRR of whitened residuals, rank by `rrpack::cv.rrr`, dense p×p whitening), `repmap_model` (RRR features→ROI patterns), `banded_ridge_model` / `banded_ridge_da_model` (encoding with `feature_sets_design`, stimulus→brain), `feature_rsa_model` (brain→features with `method = pls/pca/glmnet/ridge`, matrix targets via `feature_rsa_design$targets = F`). `pattern_model` must be positioned against these in the docs: it is the *whole-brain, spatially regularized, covariance-aware* member of the family, and it unifies the encoding and decoding directions those models treat separately. None of them compute forward patterns or spatial penalties.
- **`context$cv_spec` is populated but read by zero `fit_roi` methods.** Every shipped method pulls folds from `model$crossval`. `pattern_model` should follow the shipped convention (`crossval_samples(model$crossval, ...)`), not the documented-but-unused one.
- **Matrix-valued per-ROI payloads have exactly one exit route:** `roi_result$result`, which is dropped unless `return_fits = TRUE` (`R/mvpa_iterate.R:1327`) and is exposed in regional mode as `$fits[[i]]$predictor` (`R/regional.R:1039`). Searchlight has no equivalent. The plan keeps regional/searchlight outputs scalar via `output_schema()` (with `vector[N]` where helpful) and puts the full fit in `result$predictor` only for regional runs.
- **Surface adjacency is a shallow hole.** `neurosurf` exports `adjacency`, `neighbor_graph`, `meshToGraph`; rMVPA calls none of them. `spatial_graph()` for surfaces is a thin wrapper.
- **Two unexported volumetric adjacency builders already exist** (`build_voxel_adjacency`, `.spacenet_edges_from_feature_ids`). Promote and unify one; do not add a third.
- **Plugin test scaffolding is exported:** `mock_roi_data()`, `mock_context()`, `validate_plugin_model()`, `cv_evaluate_roi()` (`R/plugin_helpers.R`). Use them.
- **Multibasis feature identity** is already flagged as an open problem in `adocs/multibasis_prd.md:95` (duplicate voxel IDs across basis channels). The `feature_ids` + `geometry_id` contract in `spatial_graph()` resolves it for this model.
- **`.planning/architecture-refactor-prd.md` Phase 6** ("cv_labels / targets split") is the governing prior decision. `model_targets()` completes that phase rather than opening a new one.

### A.3b Positioning against existing rMVPA families

| Family | Statistical model | Direction | Spatial structure | Forward patterns | Whole-brain | What `pattern_model` reuses / adds |
|---|---|---|---|---|---|---|
| `remap_rrr_model` | RRR of whitened cross-domain residuals, `rrpack` rank CV | brain→brain (encoding→retrieval) | none | no | no (dense p×p whitening) | rank-selection idea; pattern_model adds structured noise, spatial penalties, task targets |
| `repmap_model` | RRR seed features→ROI patterns | features→brain | none | no | no | nothing directly; same "rank as diagnostic" reporting style |
| `banded_ridge_model`, `_da_model` | banded ridge encoding, per-voxel | features→brain | none | no | yes (voxelwise) | `feature_sets_design` target grouping and row weights; pattern_model adds a shared low-rank forward operator and the decoding direction from the same fit |
| `feature_rsa_model` | PLS/PCA/ridge decoding of feature vectors | brain→features | none | no | per ROI | matrix-target plumbing (`feature_rsa_design$targets`); PLS is the thresholded-PLS baseline in benchmarks |
| `spacenet_tvl1` classifier | TV-L1 linear classifier | brain→label | TV on decoding weights | Haufe post hoc | yes | edge builder and PDHG solver as 3b primitives; pattern_model regularizes forward patterns, not weights |
| **`pattern_model`** | `x = A C'y + eps`, `eps ~ (0, Psi)`, low rank, sparse + graph-smooth `A` | both, from one fit | penalties on forward patterns | model-implied (`A`), exact for calibrated scores | yes | — |

It warrants its own family because its *outputs* differ in kind: model-implied forward patterns, decoding weights, conditional-information maps, and restricted regional predictors all derived from one covariance-aware model. None of the existing families produce forward patterns or regional restriction; they are separate estimators per direction.

### A.4 Scope recommendation for v1

Build Phases 0–4 below as v1. That yields: classification and multi-response decoding, encoding predictions, sparse + Laplacian-smoothed forward patterns with diag+low-rank noise, model-implied signal-SD and conditional-information maps, restricted-ROI local performance, orthogonal rotation views, and matrix-free Haufe for the whole package. Phase 5 (confirmation inference) and Phase 6 (group) follow once v1 benchmarks are in.

---

## Part B. Implementation plan

Conventions for every phase: add each new `.R` file to `Collate` (guarded by `test_collate_field.R`); roxygen with `\code{\link{}}`, not markdown; `devtools::document()`; run the phase's tests plus `test_fit_roi.R`, `test_plugin_extension_api.R`, `test_global_analysis.R`, `test_output_schema.R`; finish each phase with `R CMD check` clean (0 errors, 0 warnings) and a fresh-context review pass. Commit per phase. Each phase should be independently mergeable.

### Phase 0 — Decisions and scaffolding (½ day)

Tasks
1. Confirm names: `pattern_model()`, class `c("pattern_model", "model_spec")`, core object class `pattern_fit`, result class `pattern_global_result`.
2. Write `R/pattern_model.R` stub with constructor + `print` + `validate_model_spec` hook so `validate_plugin_model()` passes on a trivial `fit_roi` returning baseline metrics.
3. Add a simulation helper `tests/testthat/helper-pattern_sim.R`: `sim_pattern_data(n, dims, K or q, r, snr, scenario = c("basic", "suppressor", "redundant", "diffuse", "signflip"))` returning `mvpa_dataset`, `mvpa_design`, ground-truth `A`, `C`, and the noise structure. All later known-truth tests use it.

Exit: package installs, `validate_plugin_model(pattern_model(...))` passes on the stub.

### Phase 1 — Shared foundations (1–2 days)

Independently valuable; no `pattern_model` code depends on anything beyond this phase.

1. **`model_targets()` generic** (`R/design.R`, new section) with `partition = c("train", "test")`, returning `list(values, observation_ids, response_ids, response_groups, row_weights)`.
   - `mvpa_design`: `targets` (vector or n×q matrix; fall back to `cv_labels`); new backward-compatible `targets_test` argument stored on the design, defaulting to `y_test`.
   - `feature_sets_design`: `X_train$X` / `X_test$X`, `response_groups = set`, `row_weights`.
   - `feature_rsa_design`: `F`.
   - Do **not** change `y_train.mvpa_design`.
2. **Matrix-free Haufe** (`R/importance.R`): extend `haufe_importance(W, Sigma_x = NULL, X = NULL, center = TRUE, ...)`. When `X` is supplied compute `A = X_c' Z (Z'Z)^+`, `Z = X_c W`, with optional row-block accumulation. Switch `model_importance.sda/.glmnet/.spacenet_fit` and `run_global.mvpa_model` (`global_analysis.R:321`) to the `X` path. Keep the `Sigma_x` path for callers that have one.
3. **`spatial_graph()`** (`R/pattern_spatial.R`): exported constructor over `mvpa_image_dataset` (promote `build_voxel_adjacency`; 6/18/26), `mvpa_surface_dataset` (wrap `neurosurf::adjacency`/`neighbor_graph`), `mvpa_multibasis_image_dataset` (per-channel graphs, disconnected across channels by default), and a raw `list(A=)` escape hatch. Object: `feature_ids`, sparse symmetric `A`, `degree`, `L`, `domain_type`, `geometry_id`, `n_features`. Invariant: vertex `j` ↔ column `j` of the estimator's `X`. Provide `restrict_graph(graph, keep)`.
4. **`run_global.default`** (`R/global_analysis.R`): a clear error naming `fit_roi`-style models and pointing to the class-specific method, replacing the "no applicable method" failure.

Tests: Haufe equality to the dense reference within 1e-10 on random data, including rank-deficient `W`; legacy `test_global_analysis.R` unchanged; `model_targets()` round-trips for all three designs; graph column-permutation invariance (permute features and graph together, adjacency unchanged up to relabeling); surface graph agrees with `neurosurf::adjacency` on the test mesh from `helper_surface_geom.R`.

Exit: all existing tests pass; memory of `run_global` importance path no longer scales as p².

### Phase 2 — Numerical core, no spatial penalty (2–3 days)

`R/pattern_core.R`, `R/pattern_noise.R`, `R/pattern_predict.R`.

**v1 numerical contract (fixed before coding):**

- *Objective.* `f(A, C) = 1/(2n) ||(X - Y C A') Psi^{-1/2}||_F^2 + lambda_s * sum_v ||A_v||_2 + lambda_l/2 * tr(A' L A) + lambda_2/2 ||A||_F^2`, subject to `C'C = I`. `X` and `Y` are column-centered on training rows (with `x_transform`, `y_transform`); `L` is the graph Laplacian normalized so `lambda_l` is comparable across domains (`L / mean(degree)`); `lambda_s` is expressed as a fraction of `lambda_max` (the smallest value that zeroes every row at the initial `C`), so a fixed ratio `lambda_l = rho * lambda_s` has the same meaning across folds and domains.
- *Residual covariance is fixed from a training-only pilot.* Pilot = unpenalized rank-`r_max` spectral fit on the training rows → residuals → `Psi = D + UU'` (`h` from residual spectrum). `Psi` is then held fixed for the whole alternating procedure, so the reported objective and convergence diagnostics refer to one well-defined function. An optional `noise$update = TRUE` re-estimates `Psi` once after convergence and refits (two-stage), reported as such; there is no silent per-iteration alternation.
- *Alternation.* Outer loop over (C-step, A-step) until relative objective change < tol or `max_outer`. A-step: with `C`, `Psi` fixed, the problem in `A` is convex; solved by FISTA (Phase 3) or weighted least squares (Phase 2, `lambda_s = 0`). C-step: with `A` fixed, `C* = argmin_{C'C=I} ||X~ - Y C A~'||`, `A~ = Psi^{-1/2}A`, solved by the orthogonal Procrustes closed form `C = U V'` from the SVD of `Y' X~ A~`. Orthonormal `C` removes the `A`/`C` scale ambiguity; all scale lives in `A`.
- *Initialization.* Supervised spectral: SVD of `Psi^{-1/2} X' Y (Y'Y)^{-1/2}` truncated at `r`; `C` = right singular vectors mapped back through `(Y'Y)^{-1/2}` and orthonormalized; `A` = weighted least squares given that `C`.
- *Diagnostics recorded on every fit.* Outer objective trace, inner FISTA objective trace per A-step, number of outer/inner iterations, convergence flag, `lambda_max`, `h`, pilot residual spectrum, `rank`, number of nonzero rows.

1. **Target coding** (`.pattern_encode_targets`): factor → centered one-hot (rank cap `K-1`), numeric vector/matrix → centered, optional column scaling, optional block/set weighting via a target metric. Stores `y_transform` for test-time use and the inverse for encoding output.
2. **Noise model** (`pattern_noise.R`): `Psi = D + UU'` from residuals of an initial fit; `h` by residual-spectrum heuristic (fixed default 10, user-settable, not CV-tuned in v1). Operations as closures: `apply_precision(M)`, `precision_diag()`, `restrict(keep)` (Woodbury on `D_R + U_R U_R'`), `whiten_sqrt(M)`. Shrinkage on `D` toward its median.
3. **Estimator** (`.pattern_fit(X, targets, graph = NULL, control, start = NULL)`): supervised spectral initialization (SVD of `Psi^{-1/2} X' T (T'T)^{-1/2}`), then alternate: A-step (closed-form weighted least squares in Phase 2; FISTA in Phase 3), C-step (weighted least squares then polar decomposition to enforce `C'C = I`), Psi-step every `k` outer iterations. Convergence on relative objective change. Return the `pattern_fit` object listed in the proposal (`A, C, noise, precision_A, G, target_cov, target_crosscov, x_transform, y_transform, feature_index, target_schema, rank, diagnostics`), no training images retained.
4. **Prediction** (`predict.pattern_fit(newdata, type = c("class", "prob", "decode", "encode", "scores"))`): classification via the softmax expression with priors; decoding via `Sigma_y C (I + G Phi)^{-1} u`; encoding via `A C' y`; calibrated scores `z = G^+ u` on the numerically supported subspace. Never form `G^{-1}` directly. Degenerate ROI (zero retained signal) returns priors / target means.
5. **`fit_roi.pattern_model()`** (`R/pattern_model.R`): outer folds from `crossval_samples(model$crossval, ...)` on the ROI rows; inner rank selection on inner blocked folds of the training rows by held-out decoding loss; refit selected config on outer training rows; predict outer test rows; build the **prediction ledger** (`fold_id, observation_id, prediction, truth`, matrix blocks for multi-response); metrics through `output_schema.pattern_model()` (classification: `accuracy, logloss, auc?`, `rank_selected`; regression: `r2_mean, cor_mean, rmse_mean`, `rank_selected`). Full fit rides in `result$predictor` when `return_fits = TRUE`.
6. **`run_global.pattern_model()`**: extract the domain via `get_feature_matrix()`, build/attach the graph, call the same domain evaluator once, return `pattern_global_result` with `performance_table`, `ledger`, `fold_fits` (if `return_fits`), `refit` (if `refit = TRUE`, a full-training-data descriptive fit, labeled as such), `model_spec`. Wire `mvpa_config(mode = "global")` (already routes to `run_global`).

Tests: on `sim_pattern_data("basic")` with `Psi` diagonal, recovered `A C'` matches truth in Frobenius norm within tolerance; probabilities sum to 1 and beat chance by a known margin; rank cap `K-1` enforced; encoding then decoding a noise-free `y` round-trips; `predict` on serialized/deserialized fit is identical; `cv_evaluate_roi()` from `plugin_helpers` runs the model in regional and searchlight modes; `run_global` result validates with `validate_model_spec`.

Exit: `pattern_model` works end-to-end (global, regional, searchlight) with `penalty = NULL`.

**Phase 2 as built (2026-09-07), deviations from the plan above:**

- *Target whitening.* Targets are centred **and whitened** on the training rows (`Y_w = Y_c S^{-1/2}` on the non-null eigenspace). This drops the rank deficiency of centred one-hot codes automatically (K classes give K-1 columns), makes `Yw'Yw/n = I` so the C-step is an exact orthogonal Procrustes problem, and supplies the working prior `y_w ~ N(0, I)` for decoding. `Phi = Cov(T) = I` by construction.
- *The alternation is provably at its optimum from the start.* Because the A-step's normal equation is `Psi^{-1}(A T'T - X'T)/n = 0`, `Psi` cancels and `A = X'T(T'T)^{-1}` is the exact conditional minimiser for **any** `Psi`; the C-step is exact Procrustes. Both are exact block minimisers, so the spectral initialization is already the global optimum of the unpenalized objective. Verified against 60 random restarts and an independently derived closed form (18 configurations, gaps ≤ 1.6e-13). The outer loop is retained as a convergence check and as the hook for Phase 3's penalized A-step.
- *`lambda_2` is rejected, not implemented.* Under a general `Psi` a ridge makes the A-step a Sylvester equation; the closed form `X'T(T'T + n lambda I)^{-1}` is exact only for identity noise and converges to a point that does not minimise its own objective (measured 2.5-12.4% above the true conditional minimum). `pattern_control()` therefore requires `lambda_2 = 0` until the penalized solver lands in Phase 3.
- *Initialization cost.* The reduced-rank SVD is taken on the `q_eff x p` matrix `B`, not the `n x p` fitted matrix; the right singular vectors are identical and the singular values differ by `sqrt(n)`. This avoids an `n x p` factorization at whole-brain feature counts (measured ~1500x faster at n=300, p=20000). `.estimate_pattern_noise` uses the `n x n` Gram route when `p > 2n` for the same reason.
- *Metric name.* The per-domain rank metric is `rank_mean` (the mean rank selected across folds, not necessarily an integer), not `rank_selected`.
- *Blocked CV source.* The constructor reads `design$block_var %||% design$block_var_train`, so `feature_sets_design` gets blocked CV; a design with no blocks now **warns** rather than silently using random k-fold.
- *Folds that omit a class.* Each fold's probability matrix is padded to the full class set with an exact zero column, so a class confined to one run does not crash the run; those rows score as errors with clamped log loss. Inner folds that cannot score a class are dropped from rank selection rather than propagating `NA`.
- *R² baseline.* Reported `R2` uses each fold's **training**-target mean as the baseline, matching the rank-selection loss, not the pooled test mean.
- *Not implemented in Phase 2:* target metric / feature-block weighting (`row_weights` are read from the design and a warning is issued that they are ignored), and the `graph` argument is stored but unused.

- *Repeated cross-validation is reconciled with the package convention.* Two ledgers are kept. The **fold-resolved** ledger records every prediction with its fold. The **pooled** ledger holds one record per tested observation, sorted, with repeats averaged (class probabilities averaged; continuous predictions and their training-mean baselines divided by the repeat count), exactly as `wrap_result()` in `R/mvpa_model.R` does for every other model. Metrics and the `classification_result` handed to rMVPA are built from the pooled ledger, so a row tested four times under bootstrap CV gets one vote, and the regional prediction table has unique `(roinum, .rownum)` pairs. Note the ordering is not merely a duplicate question: twofold, sequential, and bootstrap schemes are randomized, so even when each row is tested once the fold-order concatenation need not be sorted; pooling always sorts. `run_global()` returns the pooled ledger as `$ledger` and the fold-resolved one as `$fold_ledger`; ROI fits carry both as `$ledger` and `$pooled_ledger`.

**Phase 4 can rely on ledger row identity:** the pooled ledger is one sorted record per observation for every shipped CV spec.

### Phase 3 — Spatial estimator (3–4 days)

`R/pattern_spatial.R` (penalties), `R/pattern_core.R` (A-step), `R/pattern_resampling.R` (paths).

1. **A-step as FISTA**: smooth part = weighted least squares + `lambda_smooth/2 · tr(A' L A)` (matrix-free `L %*% A` with `Matrix`); prox = row-wise group soft-threshold (`lambda_sparse`) plus ridge. Step size from a power iteration on `Psi^{-1}` and `T'T/n` plus `lambda_smooth · ||L||`. Warm starts across the lambda path and across ranks.
2. **Penalty API organized by role and form.** `penalty = list(sparse = , signed_smooth = , support_smooth = )`, each either a number, `"auto"`, or `NULL`. v1 implements `sparse` (row group lasso) and `signed_smooth` (graph-Laplacian quadratic on the signed loadings, ratio `rho` to `lambda_s`, 3 ratio options). `support_smooth` (envelope `g_v` with TV) is **accepted but errors with "not implemented in v1"** so the argument's meaning is fixed now; nothing named `support_*` ever implements a quadratic penalty on signed coefficients. TV forms arrive in 3b under the same names. Manual numeric overrides bypass "auto" without enlarging the tuning grid.
3. **Column handling**: drop invalid/constant columns using training rows only, restrict the graph, keep `feature_index` so maps distinguish outside-mask / invalid / estimated-zero. No duplicate-column screening.
4. **Nested tuning** (`pattern_resampling.R`): inner blocked folds from the design's `block_var`; objective = held-out decoding loss; single parallel level (outer folds via `future`); split-scoped cache keyed on training-row identity + config.
5. **Benchmarks** (`inst/benchmarks/pattern_model/`, not run in tests). Timing covers the *whole* workflow: tuning, covariance estimation, prediction, and regional queries; never a single solver iteration. Sizes p ∈ {5k, 30k, 100k}, n = 400, r ≤ 4. Predictive comparison on identical splits with searchlight shrinkage-LDA, `spacenet_tvl1`, and thresholded PLS via `feature_rsa_model(method = "pls")`. The regional benchmark is designed to separate **a poor covariance model** from **a missed signal subspace**: on the same learned regional patterns `A_R`, compare (i) diagonal head, (ii) diag+low-rank head, (iii) training-only regional shrinkage-covariance head (an *adapted* predictor, labeled as such because it is no longer exact marginalization of the whole-brain model); then compare all three with an independently fitted, tuned shrinkage-LDA on exactly the same ROI and assessment rows. Simulations include localized nuisance modes (low-rank factors confined to a region) and redundant informative regions. Results recorded under `adocs/pattern-model-benchmarks.md`.
6. **Decision gate for C**: only if the A-step dominates at 100k voxels, add a compiled row-group prox + Laplacian product in `src/` registered in `init.c`.

Tests: objective decreases monotonically across FISTA iterations and across outer alternations; known-truth support recovery (precision/recall of nonzero rows) on `suppressor`, `redundant`, `diffuse`, `signflip` scenarios exceeds stated floors; the `signflip` scenario is a *contiguous* informative region with alternating-sign loadings, and the test records decoding loss and pattern correlation as a function of `rho` to show where signed smoothing helps and where it erases the code; the suppressor voxel has zero loading in `A` but positive conditional information (Phase 4 test reuses this); column-permutation-with-graph invariance of the fit; no assessment rows enter `x_transform`/`y_transform`/noise estimation (spy test with a poisoned test row).

Exit: v1 estimator complete, benchmarks recorded, C-kernel decision made on evidence.

**Phase 3 as built (2026-09-07), deviations from the plan above:**

- *Penalty parameterization is dimensionless.* `sparse` is a fraction of `lambda_max`, the smallest penalty that zeroes every feature at `A = 0` (an exact statement about one prox step, tested there). `signed_smooth` is the weight of the graph-Laplacian term relative to the spectral norm of the data-fit operator, after normalizing `L` to unit spectral norm. Both therefore mean the same thing across folds, ranks, and feature domains, which is what makes a fixed ratio or a short shared path legitimate. This replaces the plan's "ratio `rho` to `lambda_s`", which mixes units (the group-lasso term is linear in `A`, the Laplacian term quadratic).
- *Monotone FISTA, not plain FISTA.* The plan's test asks for a monotone inner trace, which standard FISTA does not give. The monotone variant accepts a candidate only when it does not increase the objective, so both the inner trace and the outer alternation are non-increasing by construction.
- *The inner loop does no `n x p` work.* Expanding the quadratic lets `X'T` and `tr(X Psi^{-1} X')` be computed once per outer iteration, after which each FISTA iteration costs `O(p r^2 + nnz(L) r)`. This is what makes whole-brain penalized fitting practical in R.
- *Fits are deterministic.* The power iterations that set the step size start from a fixed vector rather than the global RNG stream. Before this, two identical calls returned slightly different answers, and a leakage test could not assert exactness. The step-size estimate is inflated by 2% because power iteration converges from below and an underestimate would break the descent guarantee.
- *Memory.* A column-oriented whitening (`X W'` computed directly rather than by transposing twice) and skipping the screening copy when nothing is dropped cut peak memory at `p = 30000, n = 400` from 806 MB to 702 MB, against 92 MB for the data itself.
- *A `smooth` scenario was added to the simulator.* The plan assumed `diffuse` was the smooth counterexample to `signflip`, but `diffuse` is dense and spatially *random*, so smoothing degrades it too. A genuinely smooth blob was needed to show the contrast, and it does: smoothing raises pattern recovery from 0.86 to 0.95 on a coherent blob while lowering it from 0.99 to 0.84 on a sign-flipping one.
- *Honest finding:* on the smooth scenario, smoothing improves pattern recovery but *worsens* held-out decoding loss. Better anatomy and better prediction are not the same objective here, which is an argument for reporting both rather than tuning one and claiming the other.
- *Not implemented:* the split-scoped fit cache from the plan's item 4 (tuning is fast enough without it at the sizes measured), and parallelism inside the tuning loop (the package parallelizes at the ROI level already).

**Bugs found by review and fixed before merge (all with regression tests):**

1. *`spatial_graph()` densified.* `base::pmax` has no sparse method, so symmetrizing coerced to a dense p x p matrix: at p = 108,000 that needed 43 GB and failed. Since `signed_smooth` builds the graph automatically, whole-brain smoothing was impossible. Replaced by the sparse identity `max(a,b) = (a + b + |a - b|)/2`; the same graph now builds in 0.7 s.
2. *A missing `1/n` in the gradient operator norm.* Two consequences: `signed_smooth` was `n` times stronger than documented and scaled with the training-set size (so it was neither comparable across datasets nor stable between inner folds and the outer refit), and the A-step step size was ~n times too short. Fixed at the source; `lambda_l / data curvature` is now exactly `rho` at every `n`, and the Lipschitz estimate is within 2% of the true operator norm.
3. *Monotone FISTA treated a rejected candidate as convergence.* A rejection leaves the objective unchanged, which read as a zero relative change and stopped the solver claiming success, sometimes at the starting point. With a deliberately over-long step it returned `A = 0` and reported convergence. Now a rejection restarts the momentum, then halves the step, and never satisfies the convergence test.
4. *Convergence on the objective alone stopped too early.* The penalized objective is flat near its optimum, so the patterns, which are the scientific output, were still ~1% from the solution when the loop stopped. Convergence now requires both the objective and the iterate to settle; against a slow reference solver the fit matches to 12 digits with identical support.
5. *Multibasis ROIs collapsed onto one basis channel.* `match()` on repeated voxel ids returned channel-1 positions for every entry, so `restrict_graph` rejected the duplicate positions and *every* ROI failed; because all failed, the regional table degraded to an empty tibble with no error surfaced. Positions are now resolved by matching the k-th occurrence of an id.
6. *Tuning drew from the global RNG* when a design had no block variable, making the whole fit irreproducible on that path. Replaced by a deterministic partition.

Also: the penalty scale is now resolved per rank (so `sparse` is the same fraction of that rank's emptying penalty at every rank), the cross-rank warm start was removed (the objective is invariant under `(A, C) -> (A Q, C Q)`, so a padded lower-rank solution starts in an arbitrary rotation and was leaving the alternation unconverged), and the Laplacian's spectral norm is cached on the graph.

**Benchmarks recorded** in `adocs/pattern-model-benchmarks.md`. Headlines: a whole-brain penalized analysis with the penalty cross-validated runs in ~80 s at p = 100,000, n = 400 on one core; the restricted regional predictor matches or beats an independently fitted shrinkage LDA on the same rows in every region that carries signal; and the covariance model is not the bottleneck on this simulation (diagonal, diagonal+low-rank, and locally re-estimated heads are within 0.025 of each other). **C-kernel decision: not needed** -- the A-step is 3% of a fit at p = 100,000, because the inner loop does no `n x p` work.

**Identifiability note for Phase 4:** the objective is exactly invariant under `(A, C) -> (A Q, C Q)` for orthogonal `Q`. `A` is therefore identified only up to an `r x r` rotation, and any comparison of patterns across fits, folds, or subjects must use the column space, `A A'`, the support, or a coordinate-invariant summary such as `sqrt(diag(A Phi A'))` -- never `A` entrywise. This is the formal reason Phase 4 rotates views rather than reporting raw loadings.

### Phase 3b (optional, after benchmarks) — TV support envelope and local noise
- Row-group TV via generalized `spacenet` PDHG; support envelope `g_v` under `support_smooth`; graph-local precision term in the noise model. Same API, new `penalty` forms and `noise$type = "diag_lowrank_local"`. A regional predictor whose covariance is re-estimated locally is reported as `local_adapted`, distinct from `local_restricted` (exact marginalization of the whole-brain model).

**Provenance retained on every `pattern_fit` from Phase 2 onward** (so Phase 5 confirmation needs no estimator changes): training observation IDs, fold definition hash, `x_transform`, `y_transform`, `feature_index` + `geometry_id`, `C` and later `basis_id`, `Phi = Cov(YC)` on training rows, and the penalty/noise configuration actually used. The refit object and the fold fits are stored separately and never averaged column-wise.

### Phase 4 — Interpretation and locality (2–3 days)

`R/pattern_result.R`, `R/pattern_maps.R`.

1. **Maps** via `build_output_map()`: `model_patterns(fit, type = c("forward", "weights", "conditional_info"))` and `model_importance(fit, type = c("signal_sd", "conditional_info"))`. `signal_sd_v = sqrt([A Phi A']_vv)`; conditional information from the closed form with `h_v` from `precision_A` and `precision_diag()`. Both invariant to component-coordinate changes (tested).
2. **Empirical Haufe diagnostic**: `haufe_importance(W = fit$W, X = X_holdout)` on independent rows, returned alongside `A` with a discrepancy summary.
3. **Rotation views**: `rotate_patterns(fit, spatial = "varimax", target = "varimax")` returning `pattern_view` with `L_b, H, L_t, basis_id`, back-transforms, and a check that `L_b H L_t' == A C'`. Orthogonal only in v1.
4. **`local_performance(result, regions)`**: for each retained fold fit and region, `restrict()` the noise model, recompute `u_R, G_R`, predict the fold's test rows, score with the ledger; return a table with `whole_brain`, `local_restricted`, and (if supplied) `independent_roi` columns. Requires `return_fits = TRUE`; error with a clear message otherwise. Guarantee locality: the test-time `x_transform` must be columnwise (per-feature centering/scaling only), enforced by construction.
5. **Stability**: split-half / fold-wise principal angles between `A` column spaces → `component_stability` in the global result.
6. **Serialization**: `save_results()` support for `pattern_global_result` (maps as images, ledger and fits as RDS).

Tests: dense-Gaussian reference for regional prediction (form full `Psi_RR` explicitly on a small problem and compare to Woodbury path, 1e-10); rotation leaves `predict()` output, `signal_sd`, and conditional-info maps bit-identical; suppressor voxel: `A_v = 0` and `conditional_info_v > 0`; a regional query with a poisoned outside-region test column produces identical predictions; `local_restricted ≤ whole_brain` on average in simulation with global noise.

Exit: v1 feature-complete. Write `vignettes/Pattern_Model.Rmd` (classification + feature decoding + local access + maps) and a positioning section relative to `feature_rsa_model`, `banded_ridge_model`, `remap_rrr_model`, and searchlights.

**Phase 4 as built (2026-09-07), deviations and clarified contracts:**

- `R/pattern_maps.R` supplies original-unit forward patterns, calibrated weights, model-implied signal SD, and Gaussian conditional-information maps. Screened input columns remain `NA`. Scalar summaries are coordinate invariant; component loadings are explicitly basis dependent. Conditional information concerns continuous task scores under the working Gaussian covariance, not empirical information about categorical labels. The posterior-covariance calculation also handles singular `Phi`.
- `rotate_patterns()` returns independent orthogonal spatial/target display rotations and a coupling matrix. The base fit remains immutable; prediction and scalar-map calls delegate to it, giving bit-identical results. Display matrices use fitted feature units and whitened target coordinates; back-transforms are retained. Numerical reconstruction is checked in low-rank coordinates without allocating the full feature-by-target operator. Actual simultaneous coordinate rotations are tested separately.
- `local_performance()` accepts named regions in input matrix column positions (including logical masks). It restricts the covariance and recomputes sufficient statistics for each retained outer fit, then pools its ledger by observation before scoring. An all-screened region returns priors/target means. Independent ROI comparisons require matching fold-resolved ledgers, not unaudited scalar scores. The optional column is absent unless those ledgers are supplied. Finite-sample observed accuracy has no monotonic guarantee; the Gaussian expected-risk ordering and a simulation are tested instead.
- `pattern_haufe()` uses held-out rows and original feature units. It returns empirical and model-implied patterns plus a rotation-invariant relative Frobenius discrepancy. Rank-deficient models are compared on their identifiable score subspace. The reported score rank uses the same eigenvalue-relative cutoff as calibrated prediction; a regression reproduces the rotation-dependent matrix-entry cutoff reporting rank 2 for an effectively rank-1 decoder, and verifies the corrected rank of 1. Global results retaining fits compute fold diagnostics automatically; folds with fewer than two assessment rows are explicitly unavailable. Existing `pattern_fit` objects had no `W` field, so calibrated weights are constructed from `precision_A %*% G+` rather than assuming that field exists.
- The evaluator now retains training/assessment observation IDs, a fold hash, and a basis ID that includes the target whitening transform (identical C matrices alone do not imply the same raw target basis). IDs use `train:<design row ID>` and `test:<design row ID>` namespaces. Overlap within that namespace is rejected; independence of an external partition remains a caller obligation. This is not the complete discovery/confirmation provenance contract of Phase 5.
- `component_stability()` returns pairwise principal angles and rank-normalized projector overlap, in original feature units on the common retained columns. It reports the effective ranks and common-feature count. Empty spaces are undefined, and unequal ranks are penalized. There is no entrywise averaging or component matching. This table is attached to global results retaining at least two fits.
- `save_results.pattern_global_result()` saves the complete result as RDS and scalar refit maps as images. Multibasis image aggregation is refused by map extraction because averaging channels would change the estimand; its fitted vectors remain available and its result saves as RDS.
- `vignettes/Pattern_Model.Rmd` covers classification, continuous feature decoding, maps, regional access, rotation, diagnostics, serialization, and model-family positioning. It is linked from the pkgdown article index.
- The canonical plan and benchmark notes now live in tracked `adocs/`; historical local `.planning/` files are preserved. The full searchlight / TV-L1 / thresholded-PLS benchmark comparison remains unrun, as recorded in the benchmark note. Phase 3b and Phases 5-6 are not implemented here.

**Phase 4 validation receipt (2026-09-07):**

- Final artifact: `rMVPA_0.1.3.tar.gz`, SHA256 `10dffc175019e43e3ba0987977affb2e2c48c47760258ba8992b65b3a0693f5e`. Packaged implementation, tests, and generated help files were compared byte-for-byte with the final source manifest. `digest` is explicitly declared in Imports for provenance hashing.
- On macOS arm64, R 4.5.1: `LC_ALL=C LANG=C R_LIBS=/tmp/pattern-check-library RGL_USE_NULL=TRUE RMVPA_RUN_EXTENDED_TESTS=false R CMD check --no-manual rMVPA_0.1.3.tar.gz` finished with **0 errors, 0 warnings, 0 notes**. `carrier`, `io`, and `rrpack` plus their missing dependencies were installed in that temporary library; the user's R library was unchanged. PDF manual generation was not checked.
- The complete default package suite reported **7,266 passing assertions, zero failures, 98 test warnings, and 148 skips**. These test warnings are separate from R CMD check's WARNING count; optional/extended, CRAN-gated, and source-only checks remain skipped. The focused Phase 4 file passed **65 assertions without test warnings**. The planned integration files were also exercised, and all 32 packaged vignette HTML outputs were present; the final check rebuilt vignette outputs successfully.
- The new guide was separately rendered in a fresh R process and visually inspected, including the signal-SD/information figure. Its temporary Playwright browser was closed, ownership audit showed no automated browser processes remaining, and temporary browser tooling was removed.
- This is local evidence. Fresh-context review, hosted CI, and deployment of the new guide remain pending. The Phase 4 branch is stacked on the still-open Phase 3 PR #91; no merge or release is claimed.

### Phase 5 — Confirmation inference (3 days; separate PR)

`R/pattern_inference.R`.

1. `pattern_confirm(fit, dataset, design, block_var, inference = confirmation_plan(...))`: frozen `C` (and `basis_id`), confirmation scores `T_test = Y_test C`, unpenalized mass regression `X_test = T_test A_confirm' + E` with nuisance columns; per-voxel per-component `t`, omnibus `F`; error model options: independent, block-robust (sandwich by run), or sign-flip/permutation over blocks with max-statistic multiplicity.
2. Component tests: association (held-out correlation of `z_k` with `t_k`) and incremental value (refit reduced head, held-out loss difference, block-level resampling).
3. Provenance object: discovery/confirmation IDs with overlap check, `basis_id`, preprocessing hash, exchangeability unit.
4. `rank_supported` deferred to 5b with its own validated null.

Tests: null calibration (uniform p under label shuffling within block structure), power on `sim_pattern_data`, refusal when discovery/confirmation rows overlap.

**Phase 5 as built (2026-09-07):**

- `confirmation_plan()`, `pattern_basis()`, `pattern_confirm()`, and `pattern_component_tests()` live in `R/pattern_inference.R`. Confirmation performs unpenalized mass regression in original measurement units on frozen, training-whitened target scores with an intercept and optional numeric nuisance columns. Every input feature is tested, including discovery-screened features. Rank-deficient designs are rejected; singular sampling distributions yield unavailable tests.
- Independent Gaussian t/F tests, CR1 block sandwich approximations, and restricted-residual Rademacher block wild bootstrap are explicit choices. The bootstrap refits each component's restricted null as well as the omnibus null. It is approximate with estimated nuisance effects, not an exact permutation test. Holm families are all feature-component pairs and all feature omnibus tests separately. Bootstrap maxima span features within component (with Bonferroni across components), or all feature omnibus statistics separately. Unavailable bootstrap covariance is treated conservatively. RNG state is restored.
- Component association adjusts for nuisance but is marginal over other components. Incremental value fits full and reduced ordinary least-squares score-to-raw-target heads on caller-supplied discovery/calibration rows, then compares squared loss on untouched confirmation rows. Loss is averaged across raw target columns, within blocks, and equally across blocks. The paired t / centered wild-bootstrap inference is approximate; these are new linear heads, not the original Gaussian posterior decoder. Calibration is optional so association does not require retaining discovery measurements.
- The API requires globally meaningful discovery/confirmation IDs, explicit feature IDs, and a preprocessing recipe/unit identifier. It checks overlap and target/feature ordering, retains content and fit hashes, full frozen basis, nuisance, blocks, subject identity, and covariance. These are caller assertions about identity and independence, not proof based on legacy positional train/test prefixes. Independent errors store a small design covariance times feature residual variances; block methods store rank-by-rank covariance at each feature, never feature-by-feature covariance. Rank-1/2 block statistics are vectorized.
- The raw target basis includes training scaling and whitening, enabling exact later transport of both coefficients and covariance. Display rotations continue to use the underlying fitted inference basis. `rank_supported` remains Phase 5b: neither a component test nor selected CV rank supplies a validated sequential rank null.
- The separate prerequisite commit fixes reproduced PR #91/#93 findings: fixed-rank penalty selection evaluates the requested eligible rank, clustered ROI filtering propagates actual cluster-column positions to spatial graphs, and unavailable optional Haufe diagnostics do not discard completed evaluations. Nineteen regression assertions pass, alongside the existing integration tests. Fixed ranks are fitted directly within each inner fold, respecting its eligible rank without computing unused paths; retained test columns follow positions even when rounded centroids coincide.
- `vignettes/Pattern_Confirmation.Rmd` gives a runnable discovery/confirmation workflow, loading and component tests, calibration boundaries, covariance/provenance storage, and statistical limits. It is linked from the existing guide and pkgdown index.

**Phase 5 local evidence:**

- Independent `lm`/ANOVA and dense cluster-sandwich oracles; full coefficient/covariance rotation tests; overlap, ordering, rank deficiency, categorical rank-one and degenerate-score regressions; calibrated-head predictions and paired loss tests. The dedicated inference file passed 69 assertions without test warnings before the final package gate. Integration covers pattern files, `fit_roi`, plugin extension API, global analysis, output schema, and Collate; its sole warning is the existing two-block simulation warning.
- `inst/benchmarks/pattern_model/validate_confirmation.R`: 1,000 null experiments across 40 independently generated 30-block designs. CR1 omnibus mean p = 0.4549, rejection at nominal 5% = 8.4%; restricted wild bootstrap (199 draws) mean p = 0.4992, rejection = 4.1%. The global-null maximum family rejected in 3 of 40 designs (7.5%, a noisy estimate). This is evidence for the tested design, not universal calibration. The guide recommends the bootstrap example and explicitly reports the sandwich limitation.
- On macOS arm64, R 4.5.1 / Accelerate BLAS, the 400-row, 2,000-feature, rank-2, 40-block CR1 regression took 0.032 seconds and retained 242,584 bytes of regression summaries. The null simulation took 23.788 seconds. These are recorded workloads, not cross-platform performance guarantees.
- Full artifact package checks and rendered-guide inspection are recorded with the subsequent validation receipt; hosted CI, independent review, merging, and guide deployment remain separate gates.

**Phase 5 artifact receipt (2026-09-07):**

- The complete default `R CMD check --no-manual` on commit `0e8b1786fd7add32393f4b5e2bf731a09d0a3a8f` finished with **0 errors, 0 warnings, 0 notes**, including vignette rebuilding. The source archive SHA256 is `97afbbf363bcc231ba26fe5b8b38458d9d177128d6e0ef61183e7c1d124c8c3a`; 562 packaged code/test/help/vignette files match that commit byte-for-byte and DESCRIPTION fields match. The archive contains 33 vignette HTML outputs.
- Its default suite passed **7,346 assertions**, with zero failures, 98 test warnings, and 148 skips. These test warnings are separate from R CMD check's WARNING count; optional/extended and CRAN-gated checks remain skipped. Environment: R 4.5.1, macOS arm64, `LC_ALL=C LANG=C R_LIBS=/tmp/pattern-check-library RGL_USE_NULL=TRUE RMVPA_RUN_EXTENDED_TESTS=false`. PDF manual generation was not checked.
- Follow-up commits `7c83cfa` and `db725f9` complete fixed-rank fold handling and coincident-coordinate test-column preservation. Their owning tests and the pattern/filter suites pass. The pre-review combined Phase 5/6 head `2702c2f54ec13acbbbb9c159753809071087af90` passed **828 integration assertions**, zero failures, one existing two-block warning, and one opt-in performance skip. This explicitly distinguishes the earlier full-suite artifact from the final focused/integration evidence; the combined artifact package gate is recorded under Phase 6.
- A late review reproduced an omitted `cap_rank = TRUE` in the direct inner-fit call. The new real-fit, one-feature/two-requested-rank regression fails against `720b357` (all folds are skipped and the fallback penalty is returned). The corrected call matches independently fitted capped-fold losses; all 19 dependency regressions and the pattern/filter suites pass. The mock now enforces the real uncapped default instead of hiding it. The final post-review artifact receipt is recorded under Phase 6.
- Draft PR #94 is stacked on #93. Human review, hosted package CI, merge, and documentation deployment remain pending; no default-branch or release claim is made.

### Phase 6 — Group analysis (after Phase 5)

`R/pattern_group.R`: `pattern_group(subject_confirmations, reference_basis, spatial_mapping, effects = "random")` operating on subject loading estimates and SEs in a shared target basis; mean pattern, heterogeneity, subject expression, prediction summaries. Joint hierarchical estimator `A_s = M_s A_0 + Delta_s` is a later extension.

**Phase 6 as built (2026-09-07):**

- `pattern_group()` in `R/pattern_group.R` requires distinct subject IDs and independent confirmation observations, checks their union against every discovery set, and checks preprocessing recipe and nuisance definitions. It aligns target names and solves the exact raw-target basis relation, transporting both loadings and full coefficient covariance. Unequal target subspaces are rejected rather than Procrustes-approximated. The returned reference-coordinate hash is recomputed from the actual matrix, retaining the source discovery hash separately.
- Spatial mappings are explicit one-to-one correspondences in a common feature set, optionally selecting a shared subset. Implicit matching uses exact feature IDs. Many-to-one aggregation/interpolation is rejected because the confirmation object does not retain the cross-feature covariance required for its uncertainty. The user must perform such transformations before confirmation and establish common units, target/nuisance meaning, and conditional estimands.
- Random effects estimate the equally weighted subject population mean. Mean covariance is sample coefficient covariance divided by subject count; within-subject error is not counted twice. Component tests use t(s-1) and omnibus tests use Hotelling F(r, s-r). These are exact for iid Gaussian subject estimates with common total covariance, approximate with unequal precision or nonnormal effects. At least r+2 subjects are required. Singular distributions yield unavailable inference.
- Between-subject covariance is the PSD projection of sample covariance minus average within-subject covariance: a transparent moment estimator, not REML or a fitted joint hierarchy. Fixed effects are full-covariance GLS with asymptotic normal/chi-squared inference. The heterogeneity Q diagnostic uses r(s-1) df and treats within-subject covariance as known. It is unadjusted and distinct from the corrected component/omnibus families.
- Scalar norms, omnibus tests, heterogeneity trace, and leave-one-subject-out expression are invariant to orthogonal reference rotations. Component columns and their tests remain basis dependent. Subject expression is descriptive, not an independent group decoder. Original decoder prediction metrics remain subject-resolved and have equal-subject summaries; unavailable metrics do not silently drop a subject.
- Group covariance is transformed one feature at a time; there is no duplicated rank-by-rank-by-feature-by-subject working array. Aligned coefficients and final mean/between covariance are retained. The guide `vignettes/Pattern_Group.Rmd` covers the full public workflow, correspondence, uncertainty, heterogeneity, invariance, and the separate hierarchical/rank extensions.

**Phase 6 local evidence:**

- 54 focused group assertions cover independent GLS and Hotelling oracles, unequal-scale basis transport, orthogonal invariance, identity-safe spatial mappings, rank-one/three paths, singular covariance, covariance validity at tiny measurement scales, prediction summaries, and reference hash integrity. The two new implementation files have 96.13% line coverage from the targeted confirmation/group tests (up from 88.89% before the additional boundary and rank-three tests); this is not repository-wide coverage.
- `inst/benchmarks/pattern_model/validate_group.R`: 2,000 Gaussian sufficient-statistic null experiments with 24 subjects and rank 2 gave mean p = 0.5021 and 5.0% rejection at nominal 5%. Mean heterogeneity trace was 0.2492 versus the known 0.25. The heterogeneous-precision stress case gave mean p = 0.5031 and 4.6% rejection, supporting that tested case without making heterogeneous-precision inference exact.
- That 24-subject, 2,000-feature workload took 7.723 seconds and retained 9,224,552 bytes on the same R 4.5.1/macOS arm64 environment. Confirmation and group guides rendered and passed visual inspection, including their figures; the owned temporary browser closed and its tooling was removed.
- Phase 5 is independently committed and published as draft PR #94, stacked on #93. Phase 6 is a separate dependent change. Full package artifact receipts follow below; hosted checks, human review, merge, and guide deployment remain pending.

**Combined artifact receipt before the rank-cap review correction (2026-09-07):**

- Packaged code commit: `2702c2f54ec13acbbbb9c159753809071087af90`. Combined archive SHA256: `0b2da3f97c2c055adb4337cefe7398a469232c7267bb14872a07ad168dca1f07`. All 566 packaged code/test/help/vignette source files were compared byte-for-byte with that commit; DESCRIPTION fields also match. Later receipt-only commits change tracked planning documentation, which is excluded from the package.
- On the same R 4.5.1/macOS arm64 environment, `R CMD check --no-manual --no-tests` completed with **0 errors, 0 warnings, 0 notes**, including installation, examples, dependency/S3/help checks, and rebuilding all 34 vignette outputs. The archive reused the already-rendered Phase 5 vignette assets plus the freshly rendered group guide during `R CMD build --no-build-vignettes`; the check itself rebuilt the vignette sources successfully.
- Tests were deliberately not repeated inside this second package check: the full default core suite had already passed on the Phase 5 artifact, and the final combined head passed 828 focused/integration assertions after the two dependency follow-ups. That integration run included all 139 new confirmation, group, and review-regression assertions. It had one existing two-block warning and one opt-in performance guardrail skip. Neither a new complete default suite on the final head nor a PDF manual check is claimed.
- Draft PR #94 holds Phase 5; draft PR #95 holds Phase 6, stacked on #94. Their earlier dependencies remain #93 and #91. Hosted package CI has not run on these phase-branch PR bases (the current workflow filters master/main); human review, merge, guide deployment, and release remain pending. Local checks are not substituted for those gates.

**Final post-review receipt (2026-09-07):**

- Final packaged code commit: `b22909b3afbb217886c6d34ec172aa2aad87a387`, including the explicit capped inner-fit correction from `c999e85`. Source archive SHA256: `ae2f31696489d79a62c08b1a5eb89329752fa0c404c284a60c21908d7c8a633d`. All 566 packaged code/test/help/vignette source files and DESCRIPTION fields match this commit; the archive retains all 34 vignette HTML outputs.
- `R CMD check --no-manual --no-tests --no-vignettes` completed with **0 errors, 0 warnings, 0 notes**. Installation, examples, code, dependencies, S3 methods, help files, and vignette metadata were checked. The last one-line production correction does not alter vignette sources; all 34 had already rebuilt successfully in the preceding combined check. Vignette execution/rebuilding and the full default suite were deliberately not repeated in this final artifact check. PDF manual generation remains unchecked.
- The exact final code head separately passed **831 focused/integration assertions**, zero failures, one existing two-block warning, and one opt-in performance guardrail skip. This includes **142 new assertions**: 69 confirmation, 54 group, and 19 dependency/review regressions. The real-fit eligible-rank regression exposed the uncapped call that the original mock hid; three assertions failed before the correction and all 19 review-regression assertions pass afterward. The earlier full default suite of 7,346 assertions remains evidence for its explicitly named Phase 5 snapshot, not a claimed rerun on this head.
- Targeted inference/group line coverage remains 96.13%; those two implementation files are unchanged by the final rank-cap correction. Both guides were rendered and visually inspected, browser tooling was closed/removed, and canonical tracked plans are mirrored to the ignored local planning files. Final receipt-only commits do not change packaged code. Draft PRs #94 and #95 remain dependent on #93/#91, with human review, hosted package CI, merge, and deployment pending.

---

## Part C. File map and estimates

| File | Phase | Content |
|---|---|---|
| `R/design.R` (+) | 1 | `model_targets()` generic + 3 methods, `targets_test` |
| `R/importance.R` (+) | 1 | matrix-free `haufe_importance()` |
| `R/global_analysis.R` (+) | 1, 2 | `run_global.default`, `run_global.pattern_model` |
| `R/pattern_spatial.R` | 1, 3 | `spatial_graph()`, `restrict_graph()`, penalties, FISTA prox |
| `R/pattern_model.R` | 0, 2 | constructor, `fit_roi`, `output_schema`, `print`, preflight |
| `R/pattern_core.R` | 2, 3 | `.pattern_fit()`, init, alternating updates, `pattern_fit` object |
| `R/pattern_noise.R` | 2 | diag+low-rank noise closures, `restrict()` |
| `R/pattern_predict.R` | 2 | `predict.pattern_fit`, scores, degeneracy |
| `R/pattern_resampling.R` | 2, 3 | nested evaluation, ledger, cache |
| `R/pattern_result.R` | 2, 4 | result class, performance, `local_performance`, serialization |
| `R/pattern_maps.R` | 4 | `model_patterns`, `model_importance`, rotation views |
| `R/pattern_inference.R` | 5 | confirmation |
| `R/pattern_group.R` | 6 | group |
| `tests/testthat/helper-pattern_sim.R` | 0 | simulator |
| `tests/testthat/test_pattern_*.R` | all | one file per phase |
| `vignettes/Pattern_Model.Rmd` | 4 | user guide |
| `inst/benchmarks/pattern_model/` | 3 | benchmark scripts + recorded results |

Rough effort: Phases 0–4 ≈ 10–13 working days of implementation plus review; Phase 5 ≈ 3; Phase 6 ≈ 3–5. Phase 1 is a standalone improvement and should land first regardless.

## Part D. Open decisions for the user

1. Accept the v1 scope cut (group-lasso + Laplacian, orthogonal rotation, no rank test) with TV/envelope/local-noise as 3b?
2. Keep the numerical core inside rMVPA for now (recommended), extract to its own package only after Phase 4 stabilizes?
3. Name: `pattern_model` vs something more specific such as `rrpattern_model` / `spatial_rrr_model`. `pattern_model` is recommended for the user-facing constructor; the fit object class is `pattern_fit`.
