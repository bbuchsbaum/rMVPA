# rMVPA 0.1.2

* Added anti-leakage validator (`validate_analysis()`) with 7 cross-validation checks.
* Added permutation searchlight inference (`run_permutation_searchlight()`).
* Added feature RSA ROI connectivity outputs.
* Improved progress reporting for parallel analyses.
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
* `era_rsa_model()` now accepts combined encoding/retrieval image series split
  explicitly by `phase_var`, validates optional one-to-one item pairing, and
  supports Pearson or Spearman matched-item similarity.
* ERA-RSA can now map zero-order item-covariate correlations and adjusted,
  directional semi-partial correlations via `era_correlates`,
  `era_association`, and `era_effects`, with complete-item counts and no
  redundant unsigned R-squared maps.
* Moved 11 packages from Imports to Suggests for lazy dependency loading.
* Various bug fixes and stability improvements.
