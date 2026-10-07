# rMVPA 0.2.0: correctness, parity, speed and CRAN release plan

Status: active plan. Date: 2026-09-30. Base: `master` @ `8a13e6e`.
Scope: the release-readiness work for the whole package. The pattern_model
plan (`pattern-model-plan.md`) stays authoritative for pattern_model internals;
nothing here reopens its settled decisions.

The ordering principle is **correct → comparable → fast → clean → shipped**.
We do not optimise a path before it is pinned by characterisation tests and,
where a competitor exists, a frozen parity fixture. We do not claim to beat a
competitor without a recorded, reproducible head-to-head receipt.

---

## 0. What the audit found (2026-09-30, measured unless marked)

Environment for all measurements: Apple M2 Pro, R 4.3.3, **reference BLAS**
(`libRblas.0.dylib`), sequential. Every number is a single run with no
repetitions, so treat it as an order of magnitude only. The scripts, the profile
output and the `R CMD check` log are preserved under
`adocs/receipts/2026-09-30-audit/`.

### 0.1 Correctness defects that block any honest comparison

| ID | Finding | Evidence |
|---|---|---|
| C1 | `engine="auto"` routes **every** multiclass (≥3 levels) `mvpa_model` to SWIFT, a diagonal-Gaussian nearest-mean classifier, regardless of the requested classifier. `svmLinear` "ran" in 0.2 s with e1071 not installed. SWIFT vs legacy corclass per-centre accuracy r≈0.83. | `R/swift_searchlight.R:20-48` never checks `model$label`; `R/searchlight_engine.R:166-196` |
| C2 | Predicted class ties are broken **at random** (`max.col` default). | `R/model_fit.R:231`, `R/classifiers.R:120`, `R/performance.R:18` |
| C3 | Hyper-parameter tuning uses stratified bootstrap and **ignores run/block structure**, which leaks across runs. | `R/model_fit.R:72-77` |
| C3b | Supervised feature selection is fitted on **all outer-training labels before** inner tuning, so the inner assessment rows have already influenced the features. | `R/model_fit.R:514` (`select_features`) precedes `:544` (`tune_model`) |
| C8 | Regression `R2` is `yardstick::rsq_vec`, which is **squared correlation**, not predictive R². Offset or rescaled predictions can score as perfect. | `R/performance.R:41-55` |
| C4 | `set.seed()` called inside library code without restoring the user's RNG state. | `R/resampling_utils.R:17` (reached from `classifiers.R:558,1334`, `regional.R:370`), `spatial_nmf*.R` |
| C5 | Predicted class always = argmax of `prob()`. For SVM that means Platt probabilities (internal random CV), not decision values; e1071 also defaults `scale=TRUE`. | `R/model_fit.R:225-231`, `R/caret_models.R:107` |
| C6 | `vector_rsa` permutation p divides by `nperm+1` even when nulls are dropped as NA. | `R/vector_rsa_model.R:456-457` |
| C7 | Several perf guardrail tests compare a path with itself: the toggles they set are hard-coded `TRUE` / never read. | `mvpa_iterate.R:49`, `rsa_model.R:504`, `allgeneric.R:1050`; `test_rsa_fast_kernel_perf.R:68-73` |

### 0.2 Performance

| Scenario (100 trials, 4 classes, 5 blocked folds, r=3) | Time | ms/sphere |
|---|---|---|
| corclass, generic path, 20³ (8,000 centres) | 572 s | 71.5 |
| corclass, `auto` (→ SWIFT, *different estimator*, C1) | 1.7 s | 0.21 |
| corclass, generic path, 10³ (1,000 centres) | 67 s | 67 |
| sda_notune, generic path, 10³ | 93 s | 93 |
| dual_lda, `dual_lda_fast` (exact, C kernel), 12³ | 4.1 s | 2.4 (agent-measured) |
| rsa_model lm, generic path, 12³ | 5.9 s | 3.4 (agent-measured) |

Rprof of the generic corclass path. There are two profiles of different sizes (67 s and 205 s), and the shares vary a lot between them:

| Cost | Share of time |
|---|---|
| `futile.logger::flog.debug` evaluating namespace/`str()` *below* threshold | 25–62% |
| Per-fold `apply()` column checks (`nonzeroVarianceColumns2`) and `duplicated(t(x))` | 12–22% |
| Per-fold tibbles plus `format_result` / `merge_results` | 9–25% |
| ROI extraction (`shard_extract_roi`) | ~13% |
| `multiclass_perf` (yardstick) | ~6% |
| The actual classifier arithmetic | ~2% |

The overall conclusion is robust: overhead dominates and the arithmetic is negligible.
The ranking among the overheads is not. Phase 0's harness must re-establish it with
repeated runs before Stage A is prioritised.

For RSA, about half the time is ROI extraction, and `filter_roi` alone is about a third.

Competitor timings have **never been measured**. The only competitor scripts
(`scripts/benchmark_*_rsatoolbox.*`) need a vendored rsatoolbox that is absent.
No Python scientific stack is installed locally.

### 0.3 Parity

There are no frozen external reference outputs anywhere. The only CRAN-run parity
check compares crossnobis against a hand-written R re-implementation of
rsatoolbox (`test-crossnobis_helpers.R:88-205`).

Known definitional divergences:

| Topic | What rMVPA does |
|---|---|
| Searchlight radius | Radius is in mm with ≤ (matches nilearn; PyMVPA and CoSMoMVPA use voxels). |
| Metric aggregation | Metrics pool predictions across folds; nilearn and PyMVPA average per-fold scores. |
| AUC | Reported as `2·AUC−1`. |
| `eucdist` / `mahadist` | Not squared/P. `mahadist` estimates precision from the condition patterns themselves; rsatoolbox uses residuals. |
| RDM comparisons | No Kendall τ-a, ρ-a, cosine or whitened comparisons. |
| Searchlight null | Pooled across centres (Stelzer-style), not per-voxel. |

### 0.4 CRAN

The measured `R CMD check --as-cran` (no vignettes, R 4.3.3, GitHub neuroim2 0.17.0) gave **1 ERROR, 5 NOTEs**.

- **ERROR: 4 failing tests:**
  - an unconditional `library(spls)`;
  - the `nipalspls` / pls version;
  - a glmnet comparison at 2e-7. The cause is an API mismatch, not BLAS: the test passes `control = list(thresh, maxit)`, which glmnet 4.1.9 ignores. In a small independent reproduction, passing `thresh`/`maxit` directly cut the error from 1.6e-5 to 5.8e-9. The threshold stays.
- **Tests take about 7 min.**
- **NOTEs:**
  - `Remotes`;
  - non-CRAN Suggests `fmridesign`, `fmrilss`, `neurosurf`;
  - 23 Imports;
  - installed size 6.9 MB (a 1.7 MB benchmark CSV in `inst/extdata`);
  - diagnostic-suppressing pragmas in `src/`.
- **Clean:** R code, Rd, S3 registration, compiled-code checks.
- **External blockers:**
  - CRAN `neuroim2` is 0.13.0 while we develop against 0.17.0, and 0.17.0's exported `scale` generic breaks 3 vignettes;
  - `shard (>= 0.2.1)` but CRAN has 0.2.0.
- **Policy issues:**
  - `install_cli(dest_dir="~/.local/bin")` default;
  - 67 `\dontrun`;
  - one-sentence Description;
  - NEWS.md excluded by `.Rbuildignore`;
  - stray `tests/simple_test.R`;
  - stale `inst/CITATION`.

### 0.5 Design debt (highest-signal items)

- **Two per-ROI protocols.** One is `fit_roi`; the other is legacy `process_roi` → `internal_crossval` → `train_model`. `train_model` means "fit classifier" for `mvpa_model` but "compute whole ROI result" for RSA, MANOVA and contrast RSA.
- **Three resampling systems:**
  - modelr `resample` in `crossval.R`;
  - rsample in tuning and `create_mvpa_folds`;
  - `generate_folds`.
- **Permutation inference re-implemented 5×:**
  - `run_permutation_searchlight`;
  - vector_rsa;
  - feature_rsa;
  - banded_ridge_da;
  - spatial_nmf_inference.
- **Duplicated helpers:**
  - row correlation 3× and ridge 2×;
  - crossnobis split across two files;
  - `%||%` 4×.
- **One duplicate S3 method:** `print.rsa_design`, defined twice.
- **Too many exports.** 210 are exported, but the lifecycle registry covers 25. Exported internals include `run_future`, the `*_base` runners and `evaluate_model.*`. Exported masks: `contrasts`, `nobs`.
- **Logging on the ROOT logger.** 321 futile.logger calls use the ROOT logger, and `set_log_level` changes other packages' logging too.
- **Imports that are thin and replaceable:** stringr, memoise, tidyselect, tidyr, modelr, rsample, yardstick, and furrr/future.apply overlap.
- **Undeclared hard uses of optional packages:** c060 `epsgo` without a namespace, MGSDA, mixOmics, svd. MVPAModels' `library` field is never enforced.
- **One file too big to review:** `feature_rsa_model.R` (4,778 lines).

---

## 1. Release scope: what "core" means

Speed and parity commitments for 0.2.0 apply to the **core set** only:

| Area | Core members | Competitor reference |
|---|---|---|
| Geometry | volumetric spherical searchlight, regional (ROI) | nilearn SearchLight; PyMVPA `sphere_searchlight`; CoSMoMVPA neighbourhoods |
| Cross-validation | blocked (leave-one-run-out), k-fold with supplied indices, grouped inner tuning | sklearn LeaveOneGroupOut / GroupKFold |
| Classifiers | corclass (Pearson), nearest-centroid (Euclidean), Gaussian NB, shrinkage LDA (fixed and auto λ; `dual_lda`), linear SVM | CoSMoMVPA `cosmo_classify_*`; sklearn GaussianNB, LDA(lsqr), SVC(linear)/LinearSVC |
| Metrics | accuracy, AUC (binary, OvR macro), regression r / R² | sklearn.metrics |
| RSA | RDMs (correlation, sq. Euclidean/P, Mahalanobis/P with residual precision, crossnobis); RDM comparison (Pearson, Spearman, Kendall τ-a, cosine, regression) | rsatoolbox `calc_rdm`, `calc_rdm_crossnobis`, `compare`, `prec_from_residuals` |
| Inference | label/condition permutation for searchlight and regional, with p = (b+1)/(n+1) | sklearn `permutation_test_score`; PyMVPA MCNullDist |

Everything else (era_*, rep*, remap_rrr, naive_xdec, banded_ridge*, spatial_nmf,
feature_rsa*, contrast_rsa, pattern_model, ITEM, model-space connectivity, …)
ships as **experimental**. It must pass check and existing tests, but carries no
parity or speed claim in 0.2.0. pattern_model keeps its own plan.

---

## 2. Workstreams and phases

### Phase 0: stop the bleeding, build the yardstick

Phase 0 has no dependencies, so do it first. Each item is its own PR.

1. **Fix C1 now.** **Done** on branch `fix/searchlight-engine-estimator-identity` (commit `22116f1`):
   - `auto` never selects SWIFT. Explicit `engine="swift"` emits a message and records `attr(res, "searchlight_estimator")`.
   - New regression test file: `test_searchlight_engine_estimator_identity.R`.
   - Expected cost: multiclass models other than `dual_lda` now run on the general path (slow) until the Stage B exact engines land.
   - Of the 35 searchlight/CLI/workflow test files run, all pass except two failures that also occur on master with the same counts:
     - `test_cli.R`: 4 failures in `install_cli copies packaged wrappers`;
     - `test_mvpa_searchlight.R`: `spls` error.
   - Original option list, kept for the record:
   - SWIFT becomes eligible only for the estimator it actually computes. Register that estimator as an explicit model (e.g. `"swift_nmc"`).
   - Or SWIFT becomes eligible only when the model label maps to an *exactly equivalent* SWIFT computation, with parity proven by `expect_searchlight_parity`.
   - Add a regression test: `svmLinear` / `sda_notune` / `corclass` under `auto` must produce results identical to `engine="legacy"`.
   - Add a NEWS entry, since results from earlier versions differ.
2. **Fix C7.**
   - Delete the dead toggles, or make them real and read.
   - Every perf guardrail must compare two genuinely different code paths, or be removed.
3. **Characterisation (golden) tests for the core set.**
   - Data: small seeded synthetic data plus `inst/extdata/haxby2001_subj1`.
   - Record full outputs: per-centre metrics, per-observation predictions/probabilities, RDM vectors, p-values.
   - Store them as fixtures under `tests/testthat/fixtures/golden/`.
   - These protect every later refactor and performance change.
   - After a deliberate correctness change (Phase 1), a fixture is regenerated only with a NEWS line explaining the delta.
4. **Benchmark harness**, in `tools/bench/` (`.Rbuildignore`d, not shipped).
   - **Scenario matrix:**
     - synthetic 20³ and 40³ (≈64k centres);
     - Haxby VT (577 voxels, 96 or 864 observations);
     - one realistic whole-brain mask (~50k voxels, 200 trials, 8 runs).
     - Every core method crossed with searchlight r ∈ {mm equivalent of 2,3 voxels}, plus regional (50 ROIs).
     - Permutation with 100 permutations.
   - **Each run writes a receipt** (JSON/CSV):
     - git SHA;
     - R/BLAS/LAPACK identity and thread counts;
     - CPU;
     - wall time, ms/centre, peak RSS;
     - engine actually used;
     - output digest, so a speedup that changes results is caught.
   - **Receipts** go to `tools/bench/receipts/` and are append-only.
   - **Competitor runner:**
     - `tools/bench/competitors/` with a pinned `uv`/`pip` lock (nilearn, scikit-learn, rsatoolbox, numpy, scipy).
     - PyMVPA only if it installs cleanly on Py3; otherwise skip and say so.
     - CoSMoMVPA via Octave as optional.
     - The same data goes through NIfTI/NPZ so every tool sees identical inputs, with identical folds and thread caps.
   - **BLAS policy:**
     - Every comparison is reported at 1 thread and at N threads.
     - Every comparison is reported for R with reference BLAS and with an optimised BLAS (Accelerate/OpenBLAS), because numpy ships with OpenBLAS.
     - We claim wins only at matched BLAS and threads.
5. **Logging hot-path fix**, the smallest high-value change (measured: up to ~3× on every generic-path model).
   - Gate `flog.debug`/`flog.trace` in hot loops behind one cached logical (e.g. `.rmvpa_debug_enabled()`, refreshed when the level changes).
   - Never evaluate `paste()` / `sprintf()` arguments when disabled.
   - Move to a named `"rMVPA"` logger instead of ROOT.
   - Remove the per-batch `gc()` calls.

**Exit (release-gating part):**
- C1 and C7 fixed with tests.
- Golden fixtures in place.

**Exit (performance-track part, non-gating):**
- The harness produces repeated-run receipts, with a measured noise floor, for rMVPA and at least nilearn and rsatoolbox on the scenario matrix.
- A first scoreboard (`tools/bench/SCOREBOARD.md`) states, per scenario, rMVPA time vs best competitor. Losses are recorded as losses.

### Phase 1: correctness and parity

Phase 1 can run in parallel with Phase 3a.

**1a. Definitional fixes.** Each is a separate PR with a test and a NEWS entry.

- **C2:** deterministic tie-breaking (first maximum, matching CoSMoMVPA and sklearn argmax).
- **C3 / C3b: complete nested training.**
  - Everything learned from labels or data is fitted inside each **inner**-training split only, then refitted on the full outer-training set for the final model. This covers supervised feature selection, learned preprocessing and hyper-parameters.
  - Inner tuning uses grouped CV over `block_var` when present. Bootstrap stays available only as an explicit option.
  - Add a negative fixture: a label-dependent selector that would score above chance only if it could see inner-assessment labels must fail to do so.
- **C8: predictive R².**
  - `R2` becomes 1 − SSE/SST against a declared baseline, consistent with the pattern_model contract (fold-training mean, not the test mean).
  - Squared correlation stays available as `rsq_cor` under its own name.
  - Add fixtures with biased, rescaled and constant predictions.
- **C4:** a single internal `.with_seed()` (merge the four existing save/restore helpers). No unrestored `set.seed` anywhere in `R/`.
- **C5:**
  - An option for SVM class prediction from decision values.
  - Document `scale=TRUE`.
  - The parity fixture uses `scale=FALSE`.
- **C6: NA-aware permutation accounting**, via one shared `perm_pvalue()` helper.
  - p = (1 + #{valid nulls ≥ observed}) / (1 + #valid nulls).
  - Report the counts of requested, valid and failed permutations.
  - Return `NA`, with a recorded reason, when no valid nulls remain.
  - Optionally require a minimum valid fraction.
  - The helper standardises accounting only. Null generation stays family-specific (see Phase 3).
- **RDMs:** rsatoolbox-compatible options on the user-facing calculators:
  - squared Euclidean/P;
  - Mahalanobis/P with precision from residuals (`diag`, `shrinkage_diag`, `shrinkage_eye` to mirror `prec_from_residuals`).
  - Existing defaults stay where changing them would silently alter published analyses. The rsatoolbox-equivalent setting is documented per method.
- **RDM comparisons:** add Kendall τ-a, ρ-a, cosine and whitened (`corr_cov`/`cosine_cov`).
- **Searchlight geometry:** a `radius_units = c("mm","voxels")` argument (or a helper), so PyMVPA and CoSMoMVPA voxel-radius neighbourhoods can be reproduced. The default stays mm.
- **Metric aggregation:** a `fold_aggregation = c("pooled","mean")` option. The default stays pooled; both are documented.
- **Crossnobis:** verify unbalanced and missing-condition folds against rsatoolbox, then fix or document.

**1b. Frozen parity suite.**

- **Generator:** `data-raw/parity/` holds pinned Python scripts plus a lock file and a provenance header (versions, script hash, date). It writes small CSV/JSON files (<200 KB total) to `tests/testthat/fixtures/parity/`. CRAN tests read fixtures only. `RMVPA_PARITY_LIVE=1` regenerates them for maintainers and CI.
- **Datasets:**
  - Haxby block patterns (LORO): corclass, GaussianNB, linear SVM, fixed-λ LDA predictions, probabilities, accuracy, AUC.
  - A small anisotropic grid: neighbourhood index sets plus nilearn per-centre accuracy.
  - Haxby trials: crossnobis with residual precision.
  - Kriegeskorte92: every comparison method.
  - Deterministic sin/cos fixtures from the existing scripts, including unbalanced folds.
- **Tolerance tiers:**

  | Tier | Applies to | Requirement |
  |---|---|---|
  | Exact | fold indices, neighbourhoods, labels for corclass/dual_lda | identical |
  | Exact | RDMs, comparisons, AUC from fixed probabilities | ≤1e-10 |
  | Near | GNB, LDA posteriors, precision estimates | rtol 1e-6 |
  | Statistical | libsvm labels; sda vs sklearn auto-shrinkage | ≥98% label agreement; accuracy within 1/n |
  | Statistical | RNG-dependent nulls | compare p-value formulas exactly on fixed null vectors |

- **The "same model" rule:** a parity claim names the exact rMVPA call and the exact competitor call (estimator, parameters, folds, scoring aggregation). `vignette("parity")` publishes that table, including the documented divergences.

**Exit:** all core methods have a fixture-backed parity test that runs on CRAN in under 30 s total. The parity vignette lists every method with its tier and the settings that reproduce each toolbox.

### Phase 2: performance hill-climb

Phase 2 starts per method once that method's golden and parity tests exist.

**Rules of the climb**

- One idea per PR, with before/after receipts from the harness on the same machine. Golden and parity tests must pass unchanged, and the output digest must be equal or within the method's tier.
- Every timing is the median of ≥5 repetitions, with its IQR, after a warm-up run. The harness first measures run-to-run noise per scenario.
- Keep a change only if it wins on the scenario matrix **and** does not regress any other scenario by more than max(5%, 2× measured noise).
- The scoreboard is updated every merge. Negative results are recorded too, so ideas are not retried blindly.
- **Goal (not a release gate):** for each core scenario, rMVPA single-thread wall time ≤ 0.5× the best competitor at matched BLAS/threads, with ≤ the competitor's peak RSS. Secondary goal: generic path ≤ 3× the bare kernel cost.
- These goals are aspirations until the first competitor receipts exist. After that, the goal is re-set per scenario from the measurements.
- Each optimisation is admitted on its own receipt. The release ships whatever has merged. Performance claims in docs quote receipts only.

**Stage A: generic path scaffold.** Applies to all models, with no numerical change.

| Change | Expected effect |
|---|---|
| Logger gating (Phase 0.5) | measured ≈3× |
| Compute per-column NA/zero-variance validity **once per training split for the whole mask** (a per-voxel property, so caching it is exact), then subset per sphere. Use vectorised `colSds`/`colSums(is.na)`. | measured 1.19 → 0.08 ms/fold for the checks; ≈50% of RSA time |
| **Not output-preserving:** duplicate-column screening `duplicated(t(x))` depends on which columns share a sphere. Speed it up only in an equivalent form (e.g. hash columns once, compare within the sphere). Removing it or applying it globally is a behaviour change: a separate correctness PR with a NEWS entry. | — |
| Matrix-first ROI extraction: one T×V matrix plus integer neighbour lists, no S4 `ROIVec` per sphere (generalise the existing `.matrix_first_roi` path) | — |
| Lean results: preallocated centre×metric matrices and per-centre prediction arrays; metrics vectorised over centres; tibbles built once at the end | — |
| Fold-wise z-scoring computed once per fold, not per sphere×fold | — |

Stage A target: generic corclass ≤ 5 ms/sphere (from 67–72).

Any replacement for `filter_roi` (`R/resample.R`) must preserve:
- centre-voxel preservation;
- multibasis grouping;
- empty-ROI and min-feature refusals;
- feature identities.

The golden fixtures enforce these.

**Stage B: exact whole-brain kernels.** This is where we beat competitors.

**Eligibility contract (applies to every Stage B engine).** Each engine declares the exact regime it reproduces:
- estimator and parameters;
- preprocessing (scaling and centring inside the training fold);
- covariance or precision regime;
- priors;
- tuning (none, or a fixed grid with the inner scheme);
- outputs (labels, probabilities, decision values).

`engine="auto"` picks an engine only when the model spec falls inside that regime, and otherwise falls back to the generic path. Every engine has a registry-level eligibility test plus a parity test against the generic path on in-regime specs, and a test that out-of-regime specs fall back. A faster *different* estimator is never substituted silently; that would repeat C1.

- **Sphere-aggregation engine** (the PyMVPA GNBSearchlight idea, generalised and exact):
  - For voxel-separable sufficient statistics, one sparse centre×voxel product per fold gives every sphere at once: Σx, Σx², class-mean products, Σx·m.
  - This is **exact** corclass (Pearson), Euclidean nearest-centroid, Gaussian NB and diagonal LDA.
  - It replaces SWIFT's approximate role. SWIFT is either retired or kept as an explicitly named model.
  - Expected ~0.2–0.5 ms/sphere.
- **Native shrinkage discriminant (`sda` replacement):** maintainer priority (2026-10-01).
  - The `sda` package alone costs about 37 ms per sphere, and `sda_notune` is the recommended default.
  - Implement the same estimator in rMVPA: James–Stein shrinkage of correlations, variances and class frequencies (Schäfer–Strimmer analytic intensities), then the discriminant.
  - Golden/parity against `sda::sda` predictions and posteriors on the existing fixtures, at a stated tolerance.
  - Then run it in the incremental engine below.
- **Shrinkage LDA:**
  - Extend the `dual_lda_fast` incremental-Cholesky engine to `sda_notune`-equivalent and sklearn-equivalent shrinkage. This needs an explicit λ mapping and parity tier.
  - Route full-covariance LDA there.
  - Woodbury fold downdates instead of per-fold refits.
- **Linear SVM (low priority, maintainer 2026-10-01; after everything else in Stage B):** a compiled dual coordinate-descent solver (liblinear-style), warm-started from the neighbouring sphere along the snake order already used by dual_lda.
  - It is a **new, separately named model** (e.g. `svm_linear_l2`), with its own loss, intercept handling and multiclass strategy, matching sklearn `LinearSVC`.
  - It never accelerates the existing e1071 `svmLinear`, which is libsvm with hinge loss, an unregularised intercept and OvO.
  - Parity is against liblinear/LinearSVC.
- **RSA searchlight engine:**
  - All-centre sphere Grams in one pass give correlation, sq. Euclidean and crossnobis RDMs per centre.
  - Shared aggregation is exact only for identity or **diagonal** precision. Whole-mask whitening followed by subsetting does not reproduce ROI-specific covariance whitening. Sphere-specific full or shrinkage precision therefore stays on a per-sphere kernel (Stage C item 7). Eligibility enforces this.
  - All centres are scored with one GEMM against the precomputed model-design pseudoinverse.
  - Spearman via per-centre ranking of RDM vectors (C kernel).
- **Batched permutations:**
  - For mean-based classifiers, permuted class means are Πᵀ·X.
  - For RSA, permuting conditions permutes the model RDM, so all permutations become one GEMM.
  - `run_permutation_searchlight` stops re-running the full searchlight per permutation for eligible models.

**Stage C: compiled code** (decided by profile, not by default).

- Adopt **Rcpp + RcppArmadillo via LinkingTo** for new kernels. It adds no runtime Import burden beyond Rcpp itself. The existing C stays.
- **Candidates, in expected-value order:**
  1. fused CSR sphere-aggregation (OpenMP over centres, products computed on the fly);
  2. vectorised multiclass metrics / rank-AUC over centre blocks;
  3. per-centre Spearman ranking for RSA;
  4. liblinear-style SVM;
  5. shrinkage-LDA fit/predict via LAPACK;
  6. neighbourhood generation for the whole mask;
  7. whitened crossnobis per sphere.
- **OpenMP:** threads default to 1 under `R CMD check` (respect `OMP_THREAD_LIMIT` / `_R_CHECK_LIMIT_CORES_`), with an explicit `rMVPA.threads` option. BLAS thread control is documented to avoid oversubscription with `future` workers.

**Stage D: parallel execution.**

- The recorded CSV shows multisession slower than sequential. Re-tune chunking so each chunk carries ≥1 s of work, and prefer fork/shard where available.
- Benchmark 1/4/8 workers. Only a measured win changes a default.

**Exit:** none for the release. Phase 2 is a continuous track. At submission time:
- the scoreboard is current;
- every gap is recorded and explained;
- no golden or parity test was weakened.

### Phase 3: design consolidation (no rewrite)

**3a. Pre-CRAN, mechanical and low risk.** Can run alongside Phase 1.

- **DESCRIPTION/NAMESPACE hygiene:**
  - Declare `stats`, `utils`, `graphics`, `grDevices`.
  - Namespace or remove the classifiers that need undeclared packages (`glmnet_opt`/c060, MGSDA, mixOmics, svd).
  - Enforce MVPAModels' `library` field at load.
- **Delete:**
  - the duplicate `print.rsa_design`;
  - the empty `R/globals.R`;
  - the stub `R/run_searchlight_remap_rrr.R`;
  - the no-op `process_roi.item_model`;
  - 3 of the 4 `%||%` definitions;
  - `memo_rank`/memoise;
  - the deprecated `searchlight_mode()` and vestigial searchlight mode/profile options.

  There are no CRAN users yet, so no deprecation cycle is needed.
- **Export audit:**
  - Unexport or mark internal: `run_future`, the `*_base` runners, `evaluate_model.*`, `merge_classif_results`, `data_sample`, `sub_result`, `strip_dataset`, `crossv_*`, `prep_regional`.
  - Remove the aliases `euclidean` and `mvpa_multibasis_image_dataset`.
  - Rename `contrasts` → `contrast_matrix`, keeping a soft alias.
  - Give `nobs` the stats generic signature.
  - Extend `rmvpa_api_lifecycle()` to every exported function and every model family; experimental families get a badge in their docs.
- **Imports diet**, done only after the golden tests exist:
  - drop stringr, memoise, tidyselect, tidyr;
  - pick one of furrr / future.apply;
  - move optparse to Suggests;
  - move crayon to Suggests or plain `message()`.

  Target ≤ 16 Imports.

**3a′. Structural consolidation: desirable before 0.2.0, but not mechanical and not gating.** Each item lands only behind golden fixtures. Any item can slip past the release without blocking it.

- **Resampling unification:** one internal fold representation (index lists), used by `crossval.R`, tuning and the fast engines. This retires modelr, `create_mvpa_folds` and the rsample dependency. yardstick goes when the vectorised metrics land in Stage B, and those metrics must match the golden fixtures.
- **Shared permutation infrastructure, not one null generator.**
  - Share accounting (`perm_pvalue()`), RNG handling (`.with_seed()`), pooled vs per-voxel null aggregation and FWE options.
  - Null generation stays a family-supplied function, because exchangeability differs by family:
    - within-block label shuffles;
    - circular shifts;
    - RSA joint vs individual hypotheses;
    - fixed vs refitted representations.
  - Owners include `permutation_searchlight.R`, `vector_rsa_model.R`, `feature_rsa_model.R`, `banded_ridge_da_model.R` and `spatial_nmf_inference.R`.
  - Each family's null procedure keeps its own calibration test. Formula tests on fixed vectors do not validate a null generator.
- **One error-row constructor** for `roi_result` failures, replacing the hand-built error tibbles and the `"~"` sentinel. This is internal; output schema is unchanged.
- **Searchlight dispatch:** RSA, vector_rsa and feature_rsa fast kernels register in the engine registry (one dispatcher). `run_searchlight.default` loses its hard-coded chain.
- **Split `feature_rsa_model.R`** into spec/fit, kernels, ridge/tuning and permutation files. This is a pure move, verified by golden tests. Row-correlation and ridge helpers become shared.

**3b. Deferred past 0.2.0.** Documented so nobody starts them mid-release.

- Migrate RSA, MANOVA and contrast RSA logic into `fit_roi`, then retire `train_model` for non-classifiers and the `internal_crossval` fallback.
- Replace futile.logger with cli/rlang conditions.
- Satellite packages for spatial_nmf, banded_ridge and the shard backend. For 0.2.0 they stay in and are tagged experimental.

**Do not touch:**
- the `model_spec` / `fit_roi` / `roi_result` / `output_schema` contracts, `model$crossval` and the thin global method;
- pattern_model internals;
- ITEM's fmrilss delegation;
- the archived hrfdecoder;
- the dual_lda numerics, except behind parity tests.

### Phase 4: CRAN hardening and submission

**External, long lead time. The maintainer is assessing this (see §4, decision 1); the options below stand until that decision is made.**

- **neuroim2:** release a CRAN version that rMVPA can depend on with a version floor, and fix the exported `scale` masking.
  - **Status 2026-10-01:** neuroim2 0.20.0 is prepared (commit `b7b874c` on `feat/plot-hillclimb`; local check 0 errors, 0 warnings, 2 notes; no `scale` export) and installed locally.
  - Against it, rMVPA's golden fixtures and broad tests show no change, and the three `scale`-related vignette failures are gone.
  - After the CRAN release, set `neuroim2 (>= 0.20.0)`.
- **Vignettes using optional packages:** `Haxby_2001` fails when `randomForest`/`e1071` are absent (failed fits leave no `Accuracy` column). Every vignette must guard optional classifiers.
- **shard:** release 0.2.1 or lower the floor.
- **fmridesign, fmrilss, neurosurf:**
  - Option (a): CRAN-release them.
  - Option (b): add `Additional_repositories: https://bbuchsbaum.r-universe.dev` and guard every use with `requireNamespace`.

  Unguarded `neurosurf::` calls exist in `resample.R`, `spatial_nmf*.R` and `shard_backend.R`. Remove `Remotes`.

**In-package:**

1. **Fix the 4 test failures:**
   - `skip_if_not_installed("spls")`;
   - `pls (>= 2.9-0)` or skip by version;
   - glmnet: pass `thresh`/`maxit` directly (version-compatible), not through an ignored `control` list. Keep the 2e-7 threshold; changing any threshold needs explicit sign-off.
2. **Test budget:**
   - CRAN tier ≤ 3 min: golden + parity + fast unit tests.
   - Everything else goes behind `skip_on_cran()` / the existing extended gate.
   - CI runs the full suite on every PR. Today `.github/workflows/extended-tests.yaml` runs only `test_mvpa_searchlight.R` and `test_mvpa_regional.R`. It must be changed to run the whole extended/perf-gated set **before** more tests move behind the gate, or they stop being executed.
   - Delete `tests/simple_test.R`.
3. **Vignettes:**
   - All 34 build on CRAN dependencies within ~5 min total. Precompute or move heavy ones to pkgdown-only articles.
   - Require `albersdown (>= 2.1.0)` or guard it.
   - Fix `Haxby_2001` and the 3 neuroim2-`scale` failures.
   - Add `vignette("parity")` and `vignette("performance")`; the latter quotes scoreboard receipts with machine and BLAS.
4. **Size:**
   - Remove `inst/extdata/*_results.csv` benchmark outputs.
   - Move `inst/benchmarks/` → `tools/bench/`.
   - Installed size target < 5 MB.
5. **Policy:**
   - `install_cli` gets no default `dest_dir`.
   - Audit the 67 `\dontrun` examples individually: runnable, `\donttest`, or `\dontrun` with a stated reason (e.g. needs external data). No blanket conversion.
   - Use `message()` instead of `cat()` outside print methods.
   - Remove the diagnostic-suppressing pragmas in `src/`. The code already passes `FCONE`, and `USE_FC_LEN_T` is not currently defined. Defining it explicitly is optional, depending on the minimum supported R version (R-exts, "Fortran character strings").
   - Re-include NEWS.md in the build.
   - CITATION uses `meta$Version`.
   - Fix the dead pymvpa URL.
6. **DESCRIPTION:**
   - Rewrite Description: 3–5 sentences, single-quoted software names, DOIs (Haxby 2001; Kriegeskorte 2008; Walther 2016 for crossnobis).
   - Update Date and Version to 0.2.0.
7. **Checks**, always on the **exact submission tarball** in a clean library:
   - with only hard dependencies installed (all Suggests absent, `_R_CHECK_FORCE_SUGGESTS_=false`), which proves every optional use is conditional;
   - with all supported Suggests installed;
   - `--as-cran` on R-release and R-devel (win-builder; macOS builder; rhub Linux incl. clang-ASAN for the C/C++ code);
   - `urlchecker`;
   - `cran-comments.md`.
   - Then submit.

**Exit:** 0 ERROR / 0 WARNING; NOTEs limited to "new submission" (plus Additional_repositories if option b).

---

## 3. Sequencing and dependencies

Release gates are separate from performance research.

```
RELEASE PATH (gating)
  C1 fix + golden fixtures ─> targeted correctness fixes (C2–C8) ─┐
  release-dependency checks + 4 test failures + vignettes ────────┼─> bounded parity suite ─> Phase 4 in-package ─> submit
  Phase 3a hygiene (DESCRIPTION, deletions, exports, Imports) ────┘
  Phase 4 external (neuroim2 etc.; maintainer assessing) ─────────────────────────────────────────────────┘

PERFORMANCE TRACK (continuous, non-gating; each change admitted on its own receipt)
  harness + noise baseline ─> competitor receipts ─> Stage A ─> Stage B/C engines ─> Stage D
STRUCTURAL TRACK (non-gating): 3a′ resampling, permutation infrastructure, dispatch, file split
```

**Release critical path:**
1. C1 regression fix and golden fixtures;
2. C2–C8 correctness fixes;
3. reproducible checks against release dependencies (clean library, exact tarball, Suggests absent and present);
4. existing test failures and vignettes fixed;
5. bounded parity coverage of the core set;
6. CRAN test tier and submission.

**Start immediately, in parallel with step 1:**
- dependency compatibility;
- the four test failures;
- vignette builds;
- the extended-tests workflow fix.

The neuroim2 CRAN version is the main external risk. The performance track runs alongside: Stage A logger gating and the column-check work are cheap and likely land before release, but nothing in it blocks submission.

Each PR states the exact changes, the tests and receipts observed, any unrun checks, and the next smallest action, per `AGENTS.md`. No phase is "done" on prose. Exit criteria are checked against receipts and test output.

## 4. Maintainer decisions (2026-09-30)

1. **Non-CRAN dependencies: deferred; the maintainer is assessing them.** Status as of 2026-09-30:

   | Package | Field | CRAN status | Where rMVPA uses it |
   |---|---|---|---|
   | `fmridesign` | Suggests | not on CRAN | 1 test file only; no R/ code |
   | `fmrilss` | Suggests | not on CRAN | `R/item_design.R`, `R/item_model.R` (ITEM); 3 tests, 1 vignette |
   | `neurosurf` | Suggests | not on CRAN | surface data in `dataset.R`, `searchlight.R`, `resample.R`, `common.R`, `save_results.R`, `shard_backend.R`, `pattern_spatial.R`, `spatial_nmf*.R`; 6 tests, 1 vignette (some calls unguarded) |
   | `shard` | Suggests `(>= 0.2.1)` | CRAN has 0.2.0 (0.1.0 installed locally) | `R/shard_backend.R` (the default backend when installed); 10 tests, 2 vignettes |
   | `neuroim2` | Imports (no floor) | CRAN has 0.13.0; development uses 0.17.0 | everywhere; 0.17.0's exported `scale` generic breaks 3 vignettes |

   `albersdown`, `fmrihrf` and `fmriAR` are now on CRAN.
2. **Correctness fixes ship even though they change results** (C1–C8, with NEWS entries). RDM scaling and fold aggregation keep their current defaults, and rsatoolbox/nilearn-equivalent options are added.
3. **Rcpp + RcppArmadillo (LinkingTo) adopted** for new kernels.
4. **SWIFT will no longer be selected implicitly.** The exact sphere-aggregation engine replaces it under `engine="auto"`. SWIFT survives only as an explicitly named model if a user-facing reason for it is shown; otherwise it is retired.
5. **Experimental families stay in-package and are tagged experimental.** Each also gets a graduation track (below), so the set shrinks over time.

### 4.1 Graduation track for experimental families

A family moves from `experimental` to `stable` in `rmvpa_api_lifecycle()` when it meets all of the following:

1. it uses the `fit_roi` / `roi_result` / `output_schema` contract, with no family-specific runner, and its regional results are `regional_mvpa_result` (e.g. `run_regional.banded_ridge_model` does not meet this today);
2. it has golden characterisation fixtures, plus an independent oracle (an analytic case, a simulation with known truth, or an external implementation) at a stated tolerance;
3. it uses the shared services, not private copies: fold representation, permutation/p-value service, error-row constructor, row-correlation and ridge helpers;
4. its searchlight either has a registered engine or explicitly refuses through the registry, with no hard-coded dispatch;
5. it has a vignette that runs on CRAN dependencies within budget, with documented assumptions and limitations;
6. its tests cover refusal and edge cases (rank zero, empty ROI, single class in a fold) and finish in the CRAN tier.

Graduating a family is not on the 0.2.0 critical path. It is worked opportunistically after Phase 3a, in this order, cheapest first:

1. `contrast_rsa` / `msreve`
2. `feature_rsa`, after the file split
3. `naive_xdec`
4. `era_rsa`
5. `manova`
6. `remap_rrr`
7. `banded_ridge`, which needs runner/result-class convergence
8. the rest

The release notes list the status of every family.

---

## Appendix: review remarks (2026-09-30)

### Disposition (2026-09-30)

Every remark below was checked against source at `8a13e6e` and accepted. The plan body was revised as listed. The remarks themselves are kept verbatim.

| Remark | Verified | Where the plan changed |
|---|---|---|
| 1. Nested-training leak | `select_features` at `R/model_fit.R:514` precedes `tune_model` at `:544` | §0.1 C3b; Phase 1a C3/C3b, with a negative fixture |
| 2. glmnet API mismatch | test passes `control = list(...)` (`test_feature_rsa_ridge.R:720`); glmnet 4.1.9 installed | §0.4; Phase 4 item 1 (threshold kept) |
| 3. R² is squared correlation | `yardstick::rsq_vec` at `R/performance.R:48` | §0.1 C8; Phase 1a C8 |
| 4. Exact-engine eligibility | design review | Stage B eligibility contract; LinearSVC as a separately named model; crossnobis aggregation limited to identity or diagonal precision |
| 5. Stage A not uniformly output-preserving | `duplicated(t(train_dat))` at `R/model_fit.R:499` is per-sphere | Stage A table split; `filter_roi` invariants listed |
| 6. Permutation infrastructure | design review | Phase 1a C6 accounting spec; 3a′ shared infrastructure with family-supplied null generators |
| 7. Release gates vs research | accepted | §3 split into release path, performance track and structural track; Phase 2 target became a goal; 3a′ split from 3a |
| 8. CRAN/CI evidence | `extended-tests.yaml` runs 2 files | Phase 4 items 2 and 7 (exact tarball; Suggests absent and present) |
| Additional corrections | — | receipts preserved under `adocs/receipts/2026-09-30-audit/`; repeated-run rule; per-example `\dontrun` audit; `USE_FC_LEN_T` wording fixed |

One limit on the receipts: the second Rprof file had been overwritten by a later subagent run. Its shares differ (logger 25% vs 62%), and §0.2 now reports that range rather than a single figure.

**BEGIN REVIEW REMARKS — advisory feedback, separate from the release plan.**

These remarks record the review of the plan as presented earlier in this
conversation. They do not replace the plan or its subsequently recorded
maintainer decisions. The inspected repository was `master` at
`8a13e6ed3665b7c9fc38368cad2c33f347efe944`.

The estimator-identity rule, frozen reference fixtures, and preservation of
existing contracts are sound. The main concerns are missing correctness
requirements and a release scope that includes substantial new development.

### 1. C3 needs a complete nested-training fix

Phase 1's grouped CV change does not address supervised feature selection
happening **before** inner tuning. In [R/model_fit.R](../R/model_fit.R),
`train_model.mvpa_model()` calls `select_features()` using all outer-training
labels, then passes the selected matrix to `tune_model()`. Those labels include
inner-assessment labels.

Require feature selection and learned preprocessing to be fitted inside each
inner-training split. Add a negative fixture that detects inner-assessment
exposure.

### 2. Investigate the glmnet API mismatch before proposing tolerance changes

The proposed attribution of the glmnet failure to BLAS variance is premature.
The test in
[test_feature_rsa_ridge.R](../tests/testthat/test_feature_rsa_ridge.R) supplies
convergence settings through `control = list(thresh = 1e-14, maxit = 1e6)`,
which installed glmnet **4.1.9** does not apply.

Using the same seeded data against closed-form ridge, passing `thresh` and
`maxit` directly reduced maximum prediction error from **1.64e-5 to 5.84e-9**,
passing the existing `2e-7` comparison. This was a small independent reproduction,
not a rerun of the package test suite. Make version-compatible solver
configuration the first investigation; the result does not justify weakening
the evaluation threshold.

### 3. Add the regression R-squared discrepancy to the correctness audit

The core scope promises regression metric parity, but the current `R2` in
[R/performance.R](../R/performance.R) uses `yardstick::rsq_vec()`, which computes
squared correlation. That can report perfect performance for badly offset
predictions. Distinguish correlation from predictive R-squared, specify the
reference baseline, and add biased-prediction and constant-prediction fixtures.
Preserve the active pattern_model plan's fold-training-mean baseline contract.
[Yardstick documents the distinction between squared correlation and traditional
R-squared.](https://yardstick.tidymodels.org/reference/rsq.html)

### 4. Tighten the exact-engine eligibility contracts

Stage B needs explicit restrictions on preprocessing, covariance, priors,
tuning, and prediction outputs. In particular, the proposed liblinear-style SVM
needs its own estimator identity: it cannot silently accelerate the existing
e1071 model. LinearSVC differs in loss, intercept regularization, and multiclass
strategy. [Scikit-learn documents these
differences.](https://scikit-learn.org/stable/modules/generated/sklearn.svm.LinearSVC.html)

Likewise, shared sphere aggregation for crossnobis needs a stated precision
regime; whole-mask whitening followed by subsetting does not generally reproduce
ROI-specific covariance whitening. Require eligibility checks and a correct
fallback for configurations outside each proven regime.

### 5. Stage A is not uniformly output-preserving

The filtering optimization mixes safe caching with changed behavior. Per-column
NA/variance checks can be cached for a fixed training split, but duplicate-column
selection depends on which columns occur in the sphere. Removing that screening
changes current classifier behavior; globally applying it can discard a voxel
because of a matching voxel outside its sphere.

Preserve feature identities and classify intentional changes as correctness
changes. Also retain center-preservation, multibasis grouping, and empty-ROI
behavior when replacing `filter_roi` in [R/resample.R](../R/resample.R).

### 6. Share permutation infrastructure without imposing one null generator

Phase 3a's service must preserve family-specific exchangeability and refitting
rules. Existing code already distinguishes circular shifts, RSA joint versus
individual hypotheses, and fixed versus refitted representations. Within-block
shuffling cannot replace all of these. Relevant owners include
[R/permutation_searchlight.R](../R/permutation_searchlight.R),
[R/banded_ridge_da_model.R](../R/banded_ridge_da_model.R), and
[R/spatial_nmf_inference.R](../R/spatial_nmf_inference.R).

For C6, specify finite-null counts, behavior when no valid nulls remain, and
reporting of failed permutations. Testing the p-value formula on fixed vectors
does not validate the null-generating procedure.

### 7. Separate release gates from performance research

The reviewed critical path requires new aggregation/RSA engines and a scoreboard
target before submission. With no competitor measurements yet, "twice as fast
on every scenario with no greater memory" is an aspiration, not a justified
release gate.

Move new solvers, broad resampling consolidation, and permutation-service
migration outside the mandatory path. They also should not sit under
"mechanical and low risk." Start dependency compatibility, existing test
failures, and vignette checks immediately; admit measured optimizations
individually.

### 8. Make the CRAN and CI evidence requirements concrete

Add clean-library checks with optional dependencies absent and with supported
dependencies installed, plus checks of the exact submission tarball.
`Additional_repositories` is a valid route for optional dependencies, but
conditional use remains necessary. [CRAN policy supports that
distinction.](https://cran.r-project.org/web/packages/policies.html)

Also explicitly update the extended workflow:
[extended-tests.yaml](../.github/workflows/extended-tests.yaml) currently runs
only `test_mvpa_searchlight.R` and `test_mvpa_regional.R`. Moving additional tests
behind the extended gate could leave them unexecuted.

### Additional corrections and evidence requirements

- Preserve the audit scripts and check logs behind the "measured" claims.
- Specify repeated benchmark runs and variability before enforcing a 5%
  regression rule.
- Audit examples individually rather than blanket-converting `\dontrun`.
- The instruction to remove `USE_FC_LEN_T` is stale: it is absent from current
  `src/`, and its necessity depends on supported R versions. [R's portability
  guidance explains its
  role.](https://cran.r-project.org/doc/manuals/r-release/R-exts.html#Fortran-character-strings)

### Recommended sequence and validation status

Recommended order: **C1 and targeted correctness fixes → reproducible checks on
release dependencies → bounded parity coverage → optional measured optimizations
→ submission**.

The original review changed no files. Observed checks were source/CI inspection
and the isolated glmnet reproduction, which emitted a locale warning. Package
tests, full CRAN checks, vignettes, and competitor benchmarks were not rerun.
Appending this review changes documentation only; it does not establish a new
runtime-validation result.

The smallest recommended next action is to revise the mandatory release
checklist around these findings, then implement the C1 regression fix.

**END REVIEW REMARKS.**

---
