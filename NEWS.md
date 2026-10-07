# rMVPA (development version)

* Native SDA retains weight extraction and training-data importance for global
  and regional analyses.
* Thomaz LDA adapters accept unnamed ROI matrices and the current
  `sparsediscrim` class, probability, and score prediction interfaces.
* Debug logging uses the effective rMVPA logger threshold, including when a
  named logger has been configured.

* `install_cli()` no longer defaults `dest_dir` to `~/.local/bin`; the
  destination must be given, so nothing is written to the home directory
  unless chosen (CRAN policy).
* CRAN readiness:
  - tests requiring optional `spls`, newer `pls` (`nipalspls`) or an
    installed (non-development) rMVPA now skip when unavailable;
  - the glmnet ridge comparison passes its convergence settings as arguments
    (`glmnet` ignored the `control` list), so it now meets its original
    tolerance;
  - the `Haxby_2001` vignette omits classifiers whose optional packages are
    absent;
  - compiler-warning pragmas were removed from `src/`;
  - NEWS.md is included in the build and `tools/` excluded;
  - CITATION reads the package version;
  - a stray debug script under `tests/` was removed.
  - examples no longer use `\dontrun{}`: they run on small synthetic data,
    or sit in `\donttest{}` when slower; several previously broken examples
    (e.g. `temporal_rdm()`, `spatial_nmf_stability()`) were corrected.
* Robust `corclass` (`robust = TRUE`) computes its per-class Huber centres
  directly instead of through `neuroim2::split_reduce()`, which dispatched
  every class through `future.apply`; results are identical and fitting is
  about 4x faster.
* Hyperparameter tuning (`tune_grid`) slices the predictor matrix once per
  bootstrap split instead of rebuilding data frames for every grid row;
  splits and results are unchanged.
* `spatial_nmf_voxelwise_stats()` works with sparse volumetric maps (it called
  `neuroim2::indices()`, which has no `SparseNeuroVol` method).
* `print.manova_model()` handles one-sided formulas (`~ Y + block_var`).
* Package size: seven core vignettes ship with the package (`rMVPA`,
  `CrossValidation`, `Searchlight_Analysis`, `Regional_Analysis`, `RSA`,
  `Haxby_2001`, `Pattern_Model`); the other 27 are pkgdown articles at
  <https://bbuchsbaum.github.io/rMVPA/articles/>. Installed vignettes no
  longer embed the albersdown web fonts (~560 KB each); the website keeps
  them.
* `regression_result()` is exported, like its classification siblings.
* `gen_sample_dataset(external_test = TRUE)` no longer emits a stray message.
* **Correctness fix (changes results):** `run_searchlight()` with the default
  `engine = "auto"` no longer routes multiclass (three or more classes)
  `mvpa_model` searchlights to the SWIFT engine. SWIFT computes its own
  z-scored nearest-class-mean estimator. It was being selected regardless of
  the requested classifier, so `corclass`, `sda_notune`, `svmLinear` and other
  models silently returned SWIFT results instead of their own. `auto` now
  selects a fast engine only when it computes the specified estimator
  (currently `dual_lda_fast` for `dual_lda`). Every other classifier runs
  through the general-purpose iterator.
  - Multiclass searchlights from earlier versions run with the default engine
    should be rerun.
  - Expect longer run times for affected models until exact fast engines land.
  - SWIFT remains available through an explicit `engine = "swift"`. It then
    emits a message naming the substitution and records
    `attr(result, "searchlight_estimator") == "swift_nearest_mean"`.
* `run_regional()` for `vector_rsa_model` no longer fails when given runner
  arguments such as `verbose = FALSE` or `batch_size`. Those arguments were
  consumed and then forwarded a second time to the iterator ("formal argument
  'verbose' matched by multiple actual arguments").
* The `naive_bayes` classifier is vectorised over features: per-class
  variances use `matrixStats::colVars()`, likelihoods use one `dnorm()` call
  per class, and the softmax is row-vectorised. Results are bit-identical to
  before. A Haxby VT regional analysis went from 366 to 78 ms.
* Faster regional analysis and per-fold screening. The iterator no longer
  forces garbage collection several times per batch. A full collection walks
  the whole heap, so its cost grew with everything else in the session. Set
  `options(rMVPA.gc_each_batch = TRUE)` to restore per-batch collection. The
  per-fold zero-variance and missing-value column checks are vectorised, with
  results identical to the previous ones. On Haxby VT (leave-one-run-out), a
  `corclass` regional analysis went from 556 to 73 ms. Results are unchanged.
* Faster general-purpose searchlight and regional iteration. Debug logging is
  now checked once per run instead of on every per-fold call. Each disabled
  `futile.logger` call previously resolved the logger namespace before
  comparing thresholds, and the general path made hundreds of thousands of such
  calls. Measured on a 6x6x6 searchlight (100 trials, 5 folds, reference BLAS):
  `corclass` went from 72 to 21 ms per centre and `sda_notune` from 93 to 42 ms
  per centre. Results are unchanged.
* `set_log_level("DEBUG")` is no longer undone at the start of every run in
  interactive sessions. The internal logger setup compared a level name with a
  number, so its guard was always false.
* **Correctness fix (changes results):** predicted classes are now the exact
  maximum, with exact ties going to the first class, as in numpy and
  scikit-learn `argmax` and CoSMoMVPA. Before, `max.col()`'s default broke
  ties at random. It also treated scores within a relative 1e-5 of the maximum
  as tied, so it could return a class that was not the maximum, and results
  depended on the RNG state. This mostly affects `corclass`, whose softmax
  probabilities are nearly flat. In the Haxby VT regional fixture one of 96
  predictions changed (scissors 0.1250855 vs face 0.1250849: face had been
  chosen) and accuracy rose by one observation. A `corclass` searchlight now
  reproduces nilearn's `SearchLight` mean accuracy exactly (0.2620138889 on the
  benchmark volume) under any seed.
* New exact sphere-aggregation searchlight engine (`engine = "aggregate_fast"`),
  selected automatically for `corclass` (Pearson, mean prototypes) and
  `naive_bayes` classification searchlights. Per fold, per-voxel class means and products
  are computed once and summed over every sphere with sparse products, giving
  the same estimator as the per-sphere path: the same voxel screening,
  correlations, softmax, `zapsmall()` rounding, fold pooling and metrics.
  Centres whose aggregated values come too close to a rounding boundary are
  recomputed exactly with the per-sphere code. Accuracy and AUC maps match the
  general path at every centre (within 1e-15) on synthetic data and on Haxby
  VT, and the benchmark volume's mean accuracy matches nilearn's
  `SearchLight` exactly. The engine runs at about 0.26 ms per centre, against
  13 ms for the general path and 2.5 ms for nilearn. Data outside its regime
  (missing values, identical voxel columns, a class absent from a training
  fold) fall back to the general path automatically.
* `mvpa_model(..., class_metrics = TRUE)` works again for multiclass
  searchlight and regional analyses. The output schema did not declare the
  per-class `AUC_<class>` columns, so every ROI failed the schema width check
  and regional performance tables came back empty.
* **Metric change (changes results):** multiclass one-vs-rest AUC now ranks
  each class's own probability, as scikit-learn's
  `roc_auc_score(multi_class = "ovr")` does. The former score,
  `p_k - mean(p_-k)`, is a monotone function of `p_k` when probabilities sum
  to one, so it gives the same AUC in exact arithmetic. But rounding the mean
  turned tiny, distinct class probabilities into last-bit ties that counted
  as half credit. AUC values change only in such cases. In the Haxby VT
  regional fixtures this moved AUC by 0.2% (`sda_notune`), 0.4%
  (`naive_bayes`) and 1.6% (`dual_lda`, whose probabilities saturate). rMVPA's
  multiclass AUC now equals scikit-learn's to 15 digits on those data.
* `sda_notune`, the recommended default classifier, now fits the `sda`
  estimator natively. Shrinkage intensities, discriminant coefficients and
  posteriors match `sda::sda()` to about 1e-12 (posteriors identical after
  `sda`'s own `zapsmall()` rounding). The fit is about 10x faster and no
  longer needs the `sda`, `corpcor`, `entropy` and `fdrtool` packages. `sda`
  uses two SVDs and an eigendecomposition. The native fit uses the n x n Gram
  matrix of the centred, standardised data: the shrinkage intensity needs only
  its Frobenius norm, and the shrunk inverse correlation is applied by the
  Woodbury identity with one Cholesky factorisation (or the smaller p x p
  system when a fit has fewer voxels than training observations). A new
  searchlight engine (`engine = "sda_fast"`, selected automatically) computes
  the fit's per-voxel statistics once per fold and runs only the cross-voxel
  steps per sphere. Its maps are bit-identical to the per-sphere path's. Fits with an estimated
  correlation shrinkage of exactly zero (where `sda` uses a pseudoinverse) are
  still delegated to `sda::sda()`. Haxby VT regional analysis went from 372 to
  104 ms.
* `run_permutation_searchlight()` uses the exact searchlight engines
  (`aggregate_fast` for `corclass` and `naive_bayes`, `sda_fast` for
  `sda_notune`) under both permutation strategies with fixed folds. Data,
  neighbourhoods, folds and voxel validity are prepared once and each
  permutation only rescores. The default `"iterate"` strategy previously
  ran every permutation through the per-ROI iterator. Null distributions and
  p-values matched for the same seed in the recorded blocked-CV benchmark,
  with 15-117x faster runs for these models (6x6x6 volume). Passing an explicit
  engine other than `"auto"` bypasses preparation reuse.
* RSA searchlights (`rsa_model`) run on a new engine (`engine = "rsa_fast"`,
  selected automatically). It extracts the data once and calls
  `train_model.rsa_model()` per sphere on exactly the columns the per-ROI path
  uses, without the iterator's per-sphere ROI objects, filtering and result
  tables. Maps are identical for every distance and regression type,
  including semipartial. Spearman ranking in RDM computation now uses
  `matrixStats` (identical ranks). This also speeds up the per-ROI path.
  Haxby VT Spearman RSA searchlight (r = 6 mm): 2.13 s -> 0.23 s.
* `dual_lda` now solves problems with more features than training
  observations (most ROIs) in the dual, via the Woodbury identity: an n x n
  Cholesky instead of a p x p one. With small `gamma` the p x p system is
  severely ill-conditioned. On a 577-voxel ROI the previous solve was off by
  2e-8 relative and the dual solve by 2e-15, against a stable SVD reference.
  Predictions and metrics are unchanged on the Haxby fixtures. Haxby VT
  regional `dual_lda`: about 357 -> 126-195 ms.
* Permutation engine preparation is reused only for built-in blocked CV and
  explicit custom splits. Randomized CV keeps its existing per-permutation
  fold draws and RNG behavior. Permutation p-values now use sorted null
  lookups, preserving upper-tail ties and the +1 correction, and null results
  are concatenated once instead of copied into a growing matrix on every draw.
* RSA permutation searchlights (`rsa_model` with item permutations) also use
  the prepared-once engine path (`rsa_fast`). The permuted design's
  `item_perm` is honoured per sphere through `train_model.rsa_model()`, so
  null distributions and p-values are identical to the per-ROI path, for
  individual (correlation) and joint (`lm`) nulls under both strategies. 6-16x
  faster on a 6x6x6 volume.

# rMVPA 0.1.3

* `rsa_neural_rdm()` exposes full Pearson or Spearman neural correlation-distance
  matrices for reliability and noise-ceiling workflows. It preserves observation
  names and order, supports stimulus-mean centering, and rejects nonfinite inputs
  and constant centered patterns. Run/item exclusions remain with the caller or
  RSA design.

* Pattern-model vignettes now teach the method as a three-article path:
  `vignette("Pattern_Model")` leads with the task, the planted territories,
  and the object flow; confirmation and group articles are framed as the
  next two steps. The articles index gives them their own "Pattern Models"
  section.
* `pattern_model()` supports observation weights. Weights are read from the
  design's `row_weights` (as carried by `feature_sets_design`) or supplied
  directly via the new `weights` argument, which takes precedence. They enter
  the estimator itself -- weighted centring, weighted target whitening (the
  C-step remains an exact Procrustes problem), the residual-covariance
  estimate, and the penalized objective -- and the held-out loss that drives
  rank/penalty selection. The contract is exact: uniform weights reproduce the
  unweighted fit, integer weights are equivalent to replicating rows, and a
  zero weight is equivalent to omitting the row from training. Reported
  performance metrics remain unweighted, so weighted and unweighted runs stay
  comparable; the previous "row weights are not used" warning is gone.
* New benchmark `inst/benchmarks/pattern_model/bench_vs_baselines.R`: the
  head-to-head predictive comparison of `pattern_model` against searchlight
  shrinkage LDA (honest inner-CV sphere selection and an oracle upper bound),
  `spacenet_tvl1`, whole-brain shrinkage LDA, and CV-tuned PLS on identical
  blocked splits. Results are recorded in `adocs/pattern-model-benchmarks.md`.
* New `pattern_model()` analysis family: a pattern-first spatial reduced-rank
  model `x = A C'y + eps` with a structured residual covariance
  `Psi = D + UU'`. One fit yields condition classification, multivariate
  target decoding, encoding predictions, and forward patterns. Targets may be
  categorical or vector-valued (matrix `targets`, `feature_sets_design`, or
  `feature_rsa_design`), read through `model_targets()`. Rank is chosen by
  nested block-aware cross-validation on held-out decoding loss, with ties
  going to the smaller rank; the mean selected rank is reported as
  `rank_mean`. Predictions are kept both fold-resolved and pooled to one
  sorted record per observation with repeats averaged, matching how the rest
  of the package aggregates repeated cross-validation predictions, so metrics
  and prediction tables never double count a repeatedly tested row. Works in
  global, regional, and searchlight modes;
  `run_global()` returns a `pattern_global_result` carrying the out-of-fold
  prediction ledger, the per-fold fits, and an optional full-data refit.
  See `pattern_control()` and `predict.pattern_fit()`.
* `pattern_model()` gains spatial penalties on the forward patterns, so what
  is regularized is where task signal is expressed rather than which
  measurements help prediction. `penalty = list(sparse = )` is a row-wise
  group lasso, given as a fraction of the penalty that empties the model, so
  the same number means the same thing across folds and feature domains;
  `signed_smooth` adds graph-Laplacian smoothing of the signed loadings,
  weighted against the data-fit curvature and using a `spatial_graph()` built
  from the dataset when not supplied. Either may be `"auto"`, which
  cross-validates over a short path jointly with rank. `support_smooth` is
  reserved for a smooth support envelope and is rejected rather than quietly
  redirected, because the two encode different assumptions: smoothing signed
  loadings cannot represent a code that flips sign within a region. Penalized
  fits report `n_selected`.

* `haufe_importance()` accepts the training observations (`X = `) and computes
  activation patterns matrix-free, never forming the P x P covariance. The
  `model_importance()` methods for `sda`, `glmnet`, and `spacenet` fits and the
  averaged activation patterns of `run_global()` now use that path, so
  whole-brain global analyses no longer allocate a dense feature covariance.
* New `model_targets()` generic returns a design's model-specific targets
  (vector or matrix) with row identifiers, response identifiers, response
  groups, and row weights, for `mvpa_design`, `feature_sets_design`, and
  `feature_rsa_design`. `mvpa_design()` gains a `targets_test` argument and
  validates target row alignment. `y_train()` is unchanged.
* New `spatial_graph()` generic builds feature-aligned adjacency graphs for
  volume, multibasis, surface, and clustered datasets (plus raw adjacency
  input), with `restrict_graph()` and `graph_edges()` helpers.
* `run_global()` now reports a clear error for model classes without a global
  method instead of a generic dispatch failure.

* RSA permutation inference now states its null hypothesis. The default
  `permutation_control(rsa_null = "individual")` supports marginal
  correlations and single-predictor regressions. Regression with multiple
  design predictors, including nuisance predictors, requires explicit
  `rsa_null = "joint"`: no association with any design predictor. Raw item
  permutations destroy nuisance effects and cannot test an individual
  coefficient conditional on the others; unsupported requests now error
  before the observed searchlight runs, even if only one metric is selected.
  Results record and print the null hypothesis and identify the predictors
  covered by it.

* Circular permutation shifts now include zero independently in every block.
  This restores identity draws and partial-block shifts to the null; two-item
  blocks previously always swapped, producing a biased, potentially constant
  null pool. Single-predictor RSA diagnostics now report `NA` support for a
  constant or entirely missing RDM and warn for regression, rather than
  assigning VIF = 1 and the full item count.

* `run_permutation_searchlight()` now supports `rsa_model`. `permute_labels()`
  is an S3 generic; its `rsa_design` and `pair_rsa_design` methods relabel
  *items* within blocks rather than shuffling RDM entries, because every item
  enters `n - 1` pairs and an entry-wise shuffle would destroy that dependence
  and give an anti-conservative null. The permutation is applied to the rows of
  the neural pattern matrix in `train_model.rsa_model`, which is equivalent to
  permuting every model RDM by the inverse permutation, so the design matrix
  and the cached fast kernel are reused unchanged. `shuffle = "global"` is
  refused for designs that exclude within-block pairs. For `rsa_model` the run
  scores every model predictor against one shared null pool and returns a
  `permutation_result_set`; `metric` accepts a character vector for any model
  type, and a metric that is named but absent is now an error instead of a
  silent fallback to the first one. The per-ROI extractor also reads the
  one-row performance matrices that `rsa_model` produces, which it previously
  reduced to their first column whatever metric was requested.

* `rsa_model()` reports the effective support of its design. RDM entries share
  items, so `n_pairs` entries carry roughly `n_items` independent observations,
  and collinearity among predictor RDMs divides that further. The new
  `rsa_design_diagnostics()` computes `n_items / VIF` per predictor, the
  constructor stores it as `design_diagnostics` and prints it, and for
  `regtype = "lm"` or `"rfit"` it warns when a model predictor falls below 10
  effective items. Previously `check_collinearity` was a cliff: it stopped at
  |r| > 0.99 and said nothing at 0.95 with a dozen items, where per-ROI
  coefficients change sign through noise alone.

* The parallel runtime receipt is regenerated against current source, and the
  parallelism vignette now derives its claims from it. The receipt fingerprints
  the five files whose behaviour it measures, and that fingerprint had gone
  stale, so `test_parallelism_documentation.R` failed: the checked-in timings
  described source that no longer existed. Re-running the driver on the same
  machine, R version, and package set restores the match, with exact output
  parity across all eight scheduler/data-backend paths and identical task-frame
  sizes (13,048,800 default vs 105,712 shard bytes). Two prose claims had gone
  stale the same way and were describing a receipt two regenerations old --- one
  of them asserted shard was slower under `multisession` when the checked-in
  receipt already showed the opposite. The frame sizes, the fold reduction, the
  direction of every shard-vs-default comparison, and whether the repetition
  ranges overlap are now computed from the receipt in the vignette rather than
  written out, and the test asserts that derivation instead of asserting a
  literal number that a re-run invalidates.

* `banded_ridge_model()` scales its default alpha grid to the design instead of
  assuming one (#87). The old default, `10^seq(-2, 2, length.out = 9)`, capped
  the penalty at 100, but every solver path standardizes each training column
  and then scales band `b` by `sqrt(theta_b)`, so alpha competes with
  cross-product eigenvalues that average `(n - 1) * mean(D_b) / min(n - 1, p)`.
  That anchor grows with the number of rows and with band width: on a 1976-row,
  1500-column, five-band encoding design it is about 400, and the optimum sat
  at 1e4 — two orders of magnitude above where the grid ended. Every response
  then took the ceiling, and the pooled out-of-fold R2 was -0.300 instead of
  the -0.0075 a wider grid reached on identical data. `alphas = "auto"` (the
  new default) places nine points across six decades around that anchor, so the
  grid brackets the range in which the penalty changes the fit whatever the
  design's shape; the resolved values are on the model as `alpha_grid`. Passing
  an explicit numeric `alphas`, or a full `candidates` manifest, is unchanged.

* `run_banded_ridge()` no longer lets a truncated tuning grid pass silently
  (#87). A grid that stops below the optimum does not fail: every response takes
  the largest available alpha and the maps, metrics, and leave-one-band-out
  delta R2 come back fully populated, which is how a variance partition can be
  read off two badly fitting models. The result now carries
  `selection_diagnostics`, with one row per fitted model — full and each
  leave-one-band-out — giving the grid it could select from, the modal
  selection and its share, the shares pinned to each end, and the share
  strictly interior, plus per-model median/mean out-of-fold R2 and the share of
  responses above zero. An interior modal selection is the signature of a grid
  that contains the optimum. Two conditions now warn rather than waiting to be
  discovered. The first is a saturated grid: at least 95% of a model's
  selections taking the largest available alpha, or the smallest, or the two
  ends between them once the grid has an interior to leave empty — under the
  default per-response alpha scope a mask of signal and noise voxels sends its
  boundary mass to opposite ends, so a grid that brackets nothing can leave
  neither end near the threshold while almost nothing lands inside. The second
  is a median out-of-fold R2 below -0.05, which means the fit predicts worse
  than the mean of the data it was scored on, whatever the cause. Shares are
  computed over the alphas a response could actually have been given, so a
  fixed alpha scope is not reported as saturated, and a grid of one reports no
  boundary shares at all. `print()` on the result shows the modal alpha, the
  boundary shares, and the median R2.

* The near-zero-variance guard in `feature_rsa_model()` no longer fails ROIs
  it has no reason to fail (#85). Its threshold on the feature matrix `F` was
  absolute (`var < 1e-10`), so multiplying `F` by a constant, which leaves the
  centered PCR solution unchanged, turned a working analysis into zero usable
  ROIs; and it ran even under `feature_standardize = "center"`, where `F` is
  never divided by its column standard deviations and a constant column is
  inert. The screen now judges each column against its own magnitude, applies
  to `F` only when `"scale"` is requested, and names how many columns are
  degenerate and which is worst instead of aborting with `any()`. Judging a
  column against its own magnitude rather than against an absolute threshold
  or against the widest column in the matrix means a feature matrix that mixes
  units — the case `"scale"` exists for — is not refused, and rescaling `F`
  cannot change the verdict. Non-finite entries are reported separately and
  always. A column that is degenerate over every row now fails in
  `feature_rsa_model()` itself, where it can still be named, rather than in
  every ROI of a job whose only visible output is an all-NA performance table.
  The guard on the brain data `X`, which is always z-scored, is unchanged in
  intent but is likewise relative, so `X` in small units no longer fails.
* `feature_rsa_model()` no longer fits `pls` or `pca` components that the
  feature matrix does not define. The component count is capped at the
  numerical rank of the standardized feature matrix, for the final fit and
  separately within each segment of `ncomp_selection = "blocked"` or `"loo"`.
  Beyond that rank `pls` returns NaN or ~1e30 coefficients, which reached the
  caller as an all-NA performance table rather than as an error, and entered
  the tuning scores as if they were comparable numbers. A component count that
  some segment cannot define is no longer a candidate, so component counts are
  compared on the same segments.
* Banded-ridge tuning with the optimized solvers needs far less memory and
  skips redundant work (#84). The per-fold preprocessing receipt no longer runs
  a full `O(p^3)` direct solve just to record centering and scaling; it
  computes the same statistics directly. The solver cache verifies entries with
  a 128-bit content hash instead of retaining a copy of the standardized
  training matrix per (fold, theta) entry, which was the dominant memory
  consumer. Band Gram matrices for `dual_kernel` are cached once per training
  split under theta-free keys and shared across theta candidates, dual weight
  extraction, leave-one-band-out models, and later response chunks while the
  cache stays within its cap; the saving from Gram reuse is largest on
  reference or single-threaded BLAS builds. Inner candidates are now evaluated
  grouped by theta so each decomposition is reused across alphas within a few
  fits. Column centering, scaling, and intercept offsets on the banded-ridge
  path use recycled arithmetic instead of `sweep()`, which is bit-identical
  and avoids sweep's transposed temporaries on `n x p` matrices.
  `memory_limit_mb` now also caps the retained cache with
  least-recently-used eviction; results are unchanged at any cap. Provenance
  gains `solver_band_kernel_builds`, `solver_band_kernel_hits`,
  `solver_cache_evictions`, `solver_cache_oversize`, `solver_cache_peak_mb`,
  and `solver_cache_limit_mb`, and the work manifest gains matching per-model
  columns.
* `feature_rsa_model()` gains `feature_standardize = c("scale", "center")`.
  The feature matrix `F` has always been z-scored per training fold, which was
  undocumented and discards the variance profile of pre-reduced inputs such as
  PCA scores, making PCR component selection degenerate. `"center"` subtracts
  training-fold means only and threads through the PLS/PCA, ridge, and glmnet
  paths and their nested blocked tuning. Designs built from a similarity
  matrix `S` now default to `"center"`, because their feature matrix is
  exactly such scores; this changes results for `method = "pca"` on those
  designs (previously degenerate) and, more modestly, for PLS, ridge, and
  glmnet. Pass `feature_standardize = "scale"` to restore the old behaviour.
  For `method = "pca"` under column scaling, the constructor now warns when
  the standardized `F` has a flat correlation spectrum, and stores the
  diagnostic in `model$feature_spectrum`. The standardization contract is
  documented in a new help section (#82).
* `banded_ridge_model()` now accepts a `cross_validation` object (for example
  `blocked_cross_validation()`, `kfold_cross_validation()`, or
  `custom_cross_validation()`) for `outer_crossval`, converting it to internal
  folds with `purge` applied to the training rows, instead of misrouting it
  into the explicit-fold-list validator. `tune_crossval` is no longer
  required: the default uses at most five inner folds limited by the blocks
  (or rows, when no blocks are declared) available inside each outer-training
  set, and it is ignored when only one
  candidate survives the alpha/theta scopes, in which case inner tuning is
  skipped and reported `inner_score` values are `NA`. Inner folds and the
  scope constraints are now validated when the model is created, and the
  purge-gap error for explicit fold lists says that such lists are validated
  as supplied rather than purged (#80, #81).
* Feature RSA PLS/PCR component selection and ridge penalty selection can now
  optimize held-out `pattern_discrimination` in leakage-safe blocked inner CV;
  PLS/PCR can also optimize `pattern_rank_percentile`. MSE remains the default
  tuning objective for this release. The discrimination scorer exactly matches
  the reported correct-minus-incorrect pattern-correlation advantage while
  scaling linearly in held-out observations for a fixed voxel count and avoiding
  a quadratic candidate-correlation matrix.
* Feature RSA identification, discrimination, RDM, and permutation metrics now
  respect outer-fold candidate sets. This prevents merged out-of-fold
  predictions from being compared with targets that trained those predictions.
  Ridge can additionally tune blocked inner CV for held-out pattern rank via
  `lambda_objective = "pattern_rank_percentile"`; new effective-dimension and
  lambda-grid boundary diagnostics make over-shrinkage visible. Variance checks
  now use numerically stable two-pass kernels.
* `feature_rsa_model(method = "ridge")` adds a compact multi-response ridge
  estimator for dense feature spaces. It supports one-SVD GCV, exact analytic
  LOO, leakage-safe blocked selection, or a fixed normalized penalty; reports
  selected lambda and effective degrees of freedom rather than a component
  proxy; and has default/shard parity. Blocked tuning scores its full penalty
  path from spectral cross-products without materializing a voxel prediction
  matrix for every candidate.
* Added first-class single-domain banded-ridge encoding via
  `banded_ridge_model()` and `run_banded_ridge()`. The public workflow provides
  leakage-safe nested blocked CV, per-response or shared alpha/theta selection,
  direct/SVD/dual solvers with allocation and cache provenance, chunked spatial
  execution, exact outer-fold hyperparameter receipts, and optional OOF
  predictions or primal/dual weight retention.
* `feature_sets_design()` now stores training blocks and an explicit time-series
  declaration. Banded-ridge models reject unsafe random time-series validation,
  mismatched design rows, overlapping searchlight execution, and retained or
  intermediate allocation requests above caller-supplied limits.
* Optional `delta_sets` computes independently retuned predictive
  leave-one-band-out outer-OOF delta R2. Effects are not clipped and are not an
  additive unique/shared variance partition. The new executable
  `Banded_Ridge_Encoding` vignette includes objective scaling, output and
  storage contracts, a reproducible issue-#70 simulation audit, ecosystem
  comparison boundaries, and fixed-shape performance receipts.
* Rank-deficient RSA nuisance designs now retain aliased semi-partial terms as
  named `NA` values instead of failing or silently replacing the complete
  `sp_*` result set with missing maps.
* `era_rsa_model()` now accepts named `era_effects_block` formulas and emits
  block incremental R-squared, partial F, and rank-based numerator degrees of
  freedom from full and reduced models fit to the same complete-item set.
* Fixed a silent positive bias in `crossnobis` estimation:
  `compute_crossvalidated_means_sl()` now builds condition means from
  independent partitions instead of overlapping leave-one-partition-out
  training sets. Results produced with rMVPA 0.1.2's built-in crossnobis fold
  constructor should be recomputed; results built from independent partition
  means directly are unaffected.
* The optional shard backend now requires shard 0.2.1 or newer, which safely
  rejects invalidated shared-memory handles instead of dereferencing them.
* `run_custom_searchlight()` now gives callbacks separate training and test
  sphere matrices plus arbitrary caller-supplied `user_data`, so custom
  train/test statistics can be computed without concatenating image series.
* `mvpa_dataset()` now rejects training, test, and mask images whose spatial
  dimensions, spacing, origin, axes, or affine transforms do not agree.
* Searchlight iteration now chooses bounded, memory-aware batch sizes by
  default. Set `batch_size` explicitly to override the automatic choice.
* Custom regional and searchlight analyses now apply `.cores` for the duration
  of the call and restore the caller's previous `future` plan afterward.
* `save_results()` now writes custom searchlight metric maps as NIfTI files
  instead of storing their wrappers as auxiliary R objects.
* `save_results()` no longer loads optional surface packages merely to record
  their versions when writing a volumetric-result manifest.
* `mvpa_design()` now accepts row-aligned numeric, integer, character, logical,
  and factor vectors for blocking and splitting variables, with explicit
  alignment errors.
* Dual-LDA AUC scoring now treats differences at accumulated floating-point
  error scale as ties, preserving incremental/full searchlight parity across
  platforms without relaxing the parity threshold.
* `era_rsa_model()` now accepts combined encoding/retrieval image series split
  explicitly by `phase_var`, validates optional one-to-one item pairing, and
  supports Pearson or Spearman matched-item similarity.
* ERA-RSA can now map zero-order item-covariate correlations and adjusted,
  directional semi-partial correlations via `era_correlates`,
  `era_association`, and `era_effects`, with complete-item counts and no
  redundant unsigned R-squared maps.
* ERA-RSA item associations now use trial-specific matched-minus-nonmatch
  similarity by default, with raw matched similarity available explicitly via
  `era_association_score = "matched"`. The new `era_components` argument can
  skip identification and geometry work for association-focused searchlights.
* ERA-RSA searchlights now reuse prepared item pairing, use allocation-light
  item and RDM-vector kernels for finite data, and amortize shard dispatch with
  larger index-only batches and optional rather than unconditional batch GC.
* Eligible one-to-one volumetric ERA-RSA standard searchlights now select a
  dedicated direct-matrix engine automatically. It filters matrix columns with
  the same center-preservation rules as the general-purpose iterator, calls the same
  per-ROI scientific kernel, and supports shared-memory shard workers without
  constructing an ROI object for every sphere. The existing `engine = "legacy"`
  compatibility key requests an explicit general-iterator reference run;
  `engine = "era_rsa_fast"` requires the fast path.

# rMVPA 0.1.2

* Added anti-leakage validator (`validate_analysis()`) with 7 cross-validation checks.
* Added permutation searchlight inference (`run_permutation_searchlight()`).
* Added feature RSA ROI connectivity outputs.
* Improved progress reporting for parallel analyses.
* Moved 11 packages from Imports to Suggests for lazy dependency loading.
* Various bug fixes and stability improvements.
