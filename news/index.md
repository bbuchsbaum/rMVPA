# Changelog

## rMVPA 0.1.2

- Added anti-leakage validator
  ([`validate_analysis()`](http://bbuchsbaum.github.io/rMVPA/reference/validate_analysis.md))
  with 7 cross-validation checks.
- Added permutation searchlight inference
  ([`run_permutation_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md)).
- Added feature RSA ROI connectivity outputs.
- Improved progress reporting for parallel analyses.
- [`run_custom_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_custom_searchlight.md)
  now gives callbacks separate training and test sphere matrices plus
  arbitrary caller-supplied `user_data`, so custom train/test statistics
  can be computed without concatenating image series.
- [`mvpa_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
  now rejects training, test, and mask images whose spatial dimensions,
  spacing, origin, axes, or affine transforms do not agree.
- Searchlight iteration now chooses bounded, memory-aware batch sizes by
  default. Set `batch_size` explicitly to override the automatic choice.
- Custom regional and searchlight analyses now apply `.cores` for the
  duration of the call and restore the caller’s previous `future` plan
  afterward.
- [`save_results()`](http://bbuchsbaum.github.io/rMVPA/reference/save_results.md)
  now writes custom searchlight metric maps as NIfTI files instead of
  storing their wrappers as auxiliary R objects.
- Moved 11 packages from Imports to Suggests for lazy dependency
  loading.
- Various bug fixes and stability improvements.
