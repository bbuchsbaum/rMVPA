# Changelog

## rMVPA 0.1.2

- Added anti-leakage validator
  ([`validate_analysis()`](https://bbuchsbaum.github.io/rMVPA/reference/validate_analysis.md))
  with 7 cross-validation checks.
- Added permutation searchlight inference
  ([`run_permutation_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md)).
- Added feature RSA ROI connectivity outputs.
- Improved progress reporting for parallel analyses.
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
- Moved 11 packages from Imports to Suggests for lazy dependency
  loading.
- Various bug fixes and stability improvements.
