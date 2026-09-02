# Changelog

## rMVPA 0.1.3

- [`feature_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
  gains `feature_standardize = c("scale", "center")`. The feature matrix
  `F` has always been z-scored per training fold, which was undocumented
  and discards the variance profile of pre-reduced inputs such as PCA
  scores, making PCR component selection degenerate. `"center"`
  subtracts training-fold means only and threads through the PLS/PCA,
  ridge, and glmnet paths and their nested blocked tuning. Designs built
  from a similarity matrix `S` now default to `"center"`, because their
  feature matrix is exactly such scores; this changes results for
  `method = "pca"` on those designs (previously degenerate) and, more
  modestly, for PLS, ridge, and glmnet. Pass
  `feature_standardize = "scale"` to restore the old behaviour. For
  `method = "pca"` under column scaling, the constructor now warns when
  the standardized `F` has a flat correlation spectrum, and stores the
  diagnostic in `model$feature_spectrum`. The standardization contract
  is documented in a new help section
  ([\#82](https://github.com/bbuchsbaum/rMVPA/issues/82)).
- [`banded_ridge_model()`](https://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_model.md)
  now accepts a `cross_validation` object (for example
  [`blocked_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md),
  [`kfold_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/kfold_cross_validation.md),
  or
  [`custom_cross_validation()`](https://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md))
  for `outer_crossval`, converting it to internal folds with `purge`
  applied to the training rows, instead of misrouting it into the
  explicit-fold-list validator. `tune_crossval` is no longer required:
  the default uses at most five inner folds limited by the blocks (or
  rows, when no blocks are declared) available inside each
  outer-training set, and it is ignored when only one candidate survives
  the alpha/theta scopes, in which case inner tuning is skipped and
  reported `inner_score` values are `NA`. Inner folds and the scope
  constraints are now validated when the model is created, and the
  purge-gap error for explicit fold lists says that such lists are
  validated as supplied rather than purged
  ([\#80](https://github.com/bbuchsbaum/rMVPA/issues/80),
  [\#81](https://github.com/bbuchsbaum/rMVPA/issues/81)).
- Feature RSA PLS/PCR component selection and ridge penalty selection
  can now optimize held-out `pattern_discrimination` in leakage-safe
  blocked inner CV; PLS/PCR can also optimize `pattern_rank_percentile`.
  MSE remains the default tuning objective for this release. The
  discrimination scorer exactly matches the reported
  correct-minus-incorrect pattern-correlation advantage while scaling
  linearly in held-out observations for a fixed voxel count and avoiding
  a quadratic candidate-correlation matrix.
- Feature RSA identification, discrimination, RDM, and permutation
  metrics now respect outer-fold candidate sets. This prevents merged
  out-of-fold predictions from being compared with targets that trained
  those predictions. Ridge can additionally tune blocked inner CV for
  held-out pattern rank via
  `lambda_objective = "pattern_rank_percentile"`; new
  effective-dimension and lambda-grid boundary diagnostics make
  over-shrinkage visible. Variance checks now use numerically stable
  two-pass kernels.
- `feature_rsa_model(method = "ridge")` adds a compact multi-response
  ridge estimator for dense feature spaces. It supports one-SVD GCV,
  exact analytic LOO, leakage-safe blocked selection, or a fixed
  normalized penalty; reports selected lambda and effective degrees of
  freedom rather than a component proxy; and has default/shard parity.
  Blocked tuning scores its full penalty path from spectral
  cross-products without materializing a voxel prediction matrix for
  every candidate.
- Added first-class single-domain banded-ridge encoding via
  [`banded_ridge_model()`](https://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_model.md)
  and
  [`run_banded_ridge()`](https://bbuchsbaum.github.io/rMVPA/reference/run_banded_ridge.md).
  The public workflow provides leakage-safe nested blocked CV,
  per-response or shared alpha/theta selection, direct/SVD/dual solvers
  with allocation and cache provenance, chunked spatial execution, exact
  outer-fold hyperparameter receipts, and optional OOF predictions or
  primal/dual weight retention.
- [`feature_sets_design()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md)
  now stores training blocks and an explicit time-series declaration.
  Banded-ridge models reject unsafe random time-series validation,
  mismatched design rows, overlapping searchlight execution, and
  retained or intermediate allocation requests above caller-supplied
  limits.
- Optional `delta_sets` computes independently retuned predictive
  leave-one-band-out outer-OOF delta R2. Effects are not clipped and are
  not an additive unique/shared variance partition. The new executable
  `Banded_Ridge_Encoding` vignette includes objective scaling, output
  and storage contracts, a reproducible issue-#70 simulation audit,
  ecosystem comparison boundaries, and fixed-shape performance receipts.
- Rank-deficient RSA nuisance designs now retain aliased semi-partial
  terms as named `NA` values instead of failing or silently replacing
  the complete `sp_*` result set with missing maps.
- [`era_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md)
  now accepts named `era_effects_block` formulas and emits block
  incremental R-squared, partial F, and rank-based numerator degrees of
  freedom from full and reduced models fit to the same complete-item
  set.
- Fixed a silent positive bias in `crossnobis` estimation:
  [`compute_crossvalidated_means_sl()`](https://bbuchsbaum.github.io/rMVPA/reference/compute_crossvalidated_means_sl.md)
  now builds condition means from independent partitions instead of
  overlapping leave-one-partition-out training sets. Results produced
  with rMVPA 0.1.2’s built-in crossnobis fold constructor should be
  recomputed; results built from independent partition means directly
  are unaffected.
- The optional shard backend now requires shard 0.2.1 or newer, which
  safely rejects invalidated shared-memory handles instead of
  dereferencing them.
- [`run_custom_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_custom_searchlight.md)
  now gives callbacks separate training and test sphere matrices plus
  arbitrary caller-supplied `user_data`, so custom train/test statistics
  can be computed without concatenating image series.
- [`mvpa_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
  now rejects training, test, and mask images whose spatial dimensions,
  spacing, origin, axes, or affine transforms do not agree.
- Searchlight iteration now chooses bounded, memory-aware batch sizes by
  default. Set `batch_size` explicitly to override the automatic choice.
- Custom regional and searchlight analyses now apply `.cores` for the
  duration of the call and restore the caller’s previous `future` plan
  afterward.
- [`save_results()`](https://bbuchsbaum.github.io/rMVPA/reference/save_results.md)
  now writes custom searchlight metric maps as NIfTI files instead of
  storing their wrappers as auxiliary R objects.
- [`save_results()`](https://bbuchsbaum.github.io/rMVPA/reference/save_results.md)
  no longer loads optional surface packages merely to record their
  versions when writing a volumetric-result manifest.
- [`mvpa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
  now accepts row-aligned numeric, integer, character, logical, and
  factor vectors for blocking and splitting variables, with explicit
  alignment errors.
- Dual-LDA AUC scoring now treats differences at accumulated
  floating-point error scale as ties, preserving incremental/full
  searchlight parity across platforms without relaxing the parity
  threshold.
- [`era_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md)
  now accepts combined encoding/retrieval image series split explicitly
  by `phase_var`, validates optional one-to-one item pairing, and
  supports Pearson or Spearman matched-item similarity.
- ERA-RSA can now map zero-order item-covariate correlations and
  adjusted, directional semi-partial correlations via `era_correlates`,
  `era_association`, and `era_effects`, with complete-item counts and no
  redundant unsigned R-squared maps.
- ERA-RSA item associations now use trial-specific
  matched-minus-nonmatch similarity by default, with raw matched
  similarity available explicitly via
  `era_association_score = "matched"`. The new `era_components` argument
  can skip identification and geometry work for association-focused
  searchlights.
- ERA-RSA searchlights now reuse prepared item pairing, use
  allocation-light item and RDM-vector kernels for finite data, and
  amortize shard dispatch with larger index-only batches and optional
  rather than unconditional batch GC.
- Eligible one-to-one volumetric ERA-RSA standard searchlights now
  select a dedicated direct-matrix engine automatically. It filters
  matrix columns with the same center-preservation rules as the
  general-purpose iterator, calls the same per-ROI scientific kernel,
  and supports shared-memory shard workers without constructing an ROI
  object for every sphere. The existing `engine = "legacy"`
  compatibility key requests an explicit general-iterator reference run;
  `engine = "era_rsa_fast"` requires the fast path.

## rMVPA 0.1.2

- Added anti-leakage validator
  ([`validate_analysis()`](https://bbuchsbaum.github.io/rMVPA/reference/validate_analysis.md))
  with 7 cross-validation checks.
- Added permutation searchlight inference
  ([`run_permutation_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md)).
- Added feature RSA ROI connectivity outputs.
- Improved progress reporting for parallel analyses.
- Moved 11 packages from Imports to Suggests for lazy dependency
  loading.
- Various bug fixes and stability improvements.
