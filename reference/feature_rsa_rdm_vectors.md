# Extract Per-ROI Predicted and Observed RDM Vectors from Feature RSA Results

Convenience helper to pull the compact lower-triangle predicted (and
optionally observed) RDM vectors stored by
`feature_rsa_model(..., return_rdm_vectors = TRUE)` from a
`regional_mvpa_result`. Results can come either from in-memory `$fits`
or from file-backed batches written via
`run_regional(..., save_rdm_vectors_dir = ...)`.

## Usage

``` r
feature_rsa_rdm_vectors(x)
```

## Arguments

- x:

  A `regional_mvpa_result` returned by
  [`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
  for a `feature_rsa_model`, or a tibble/data frame with columns
  `roinum` and `rdm_vec`. A `regional_mvpa_result` may store vectors
  either in-memory or on disk in `$rdm_batch_dir`.

## Value

A tibble with one row per ROI and columns:

- roinum:

  ROI id.

- n_obs:

  Number of observations contributing to the vector.

- observation_index:

  List-column of observation ordering used for the predicted RDM.

- rdm_vec:

  List-column containing the lower-triangle predicted RDM vector for
  that ROI.

- observed_rdm_vec:

  List-column containing the lower-triangle observed RDM vector for that
  ROI (if available).

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_regional(
  feature_rsa_model(dataset, design, method = "pls", return_rdm_vectors = TRUE),
  region_mask
)
vecs <- feature_rsa_rdm_vectors(res)
} # }
```
