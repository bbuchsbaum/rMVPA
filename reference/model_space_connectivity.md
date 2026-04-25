# Model-space representational connectivity from a fitted rMVPA result

Unified entry point for second-order representational connectivity.
Accepts a matrix/list/tibble of per-unit RDM vectors, or a fitted
regional result carrying per-ROI model-space fingerprints / RDM vectors,
and forwards to
[`rdm_model_space_connectivity`](http://bbuchsbaum.github.io/rMVPA/reference/rdm_model_space_connectivity.md)
for the projection and similarity summaries.

## Usage

``` r
model_space_connectivity(result, model_rdms = NULL, ...)

# Default S3 method
model_space_connectivity(
  result,
  model_rdms = NULL,
  method = c("pearson", "spearman"),
  basis = c("pca", "qr"),
  use = c("complete.obs", "everything"),
  tol = 1e-08,
  return_projected = FALSE,
  ...
)

# S3 method for class 'regional_mvpa_result'
model_space_connectivity(
  result,
  model_rdms = NULL,
  method = c("pearson", "spearman"),
  basis = c("pca", "qr"),
  use = c("complete.obs", "everything"),
  tol = 1e-08,
  return_projected = FALSE,
  prefer = c("fingerprint", "rdm"),
  ...
)

# S3 method for class 'searchlight_result'
model_space_connectivity(
  result,
  model_rdms = NULL,
  k = 20L,
  scale = c("norm", "raw"),
  seeds = c("kmeans"),
  nstart = 10L,
  iter.max = 30L,
  random_seed = NULL,
  build_maps = TRUE,
  ...
)
```

## Arguments

- result:

  Either a numeric matrix / list / tibble of per-unit RDM vectors
  (forwarded directly to `rdm_model_space_connectivity`), or a fitted
  regional rMVPA object that carries fingerprints / RDM vectors per ROI
  (e.g. `regional_mvpa_result` from
  [`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
  on an `rsa_model(..., return_fingerprint = TRUE)` or a
  `feature_rsa_model(..., return_rdm_vectors = TRUE)`).

- model_rdms:

  Named list of square symmetric model RDMs / `dist` objects, or a
  numeric matrix with RDM cells in rows. Required when the per-unit
  summaries are RDM vectors. Optional when fingerprints already live on
  `result`: in that case the fingerprints are taken as the ROI-by-axis
  matrix `F` and the function returns `tcrossprod(F)` plus its
  decompositions, without re-projecting.

- ...:

  Additional arguments forwarded to underlying methods.

- method, basis, use, tol, return_projected:

  Forwarded to
  [`rdm_model_space_connectivity()`](http://bbuchsbaum.github.io/rMVPA/reference/rdm_model_space_connectivity.md).
  See that function for details.

- prefer:

  Either `"fingerprint"` (default) or `"rdm"` when both representations
  are available on `result`.

- k:

  Integer; number of anchors to select. Default `20`, clamped to the
  number of available searchlight centers.

- scale:

  One of `"norm"` (default; row-normalize fingerprints so similarity is
  cosine) or `"raw"` (use raw fingerprint inner products).

- seeds:

  Anchor selection strategy: `"kmeans"` (default; cluster fingerprints
  and pick the searchlight closest to each centroid) or an integer
  vector of explicit center IDs to use as anchors.

- nstart, iter.max:

  Forwarded to [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).
  Default `nstart = 10`, `iter.max = 30`.

- random_seed:

  Integer seed for reproducible k-means starts. Default `NULL` (no
  seeding). When supplied, the caller's global RNG state is restored
  after anchor selection.

- build_maps:

  Logical; if `TRUE` (default) build one `NeuroVol`/`NeuroSurface` per
  anchor via
  [`build_output_map`](http://bbuchsbaum.github.io/rMVPA/reference/build_output_map.md).
  Skip with `FALSE` if you only need the numerical similarity matrix.

## Value

An object of class `rdm_model_space_connectivity` (see
[`rdm_model_space_connectivity()`](http://bbuchsbaum.github.io/rMVPA/reference/rdm_model_space_connectivity.md)
for details).

For `searchlight_result` inputs, an object of class
`model_space_anchor_connectivity`. This is an anchor summary, not a full
searchlight-by-searchlight matrix. It contains an `n_centers x k`
`similarity` matrix, anchor center IDs, optional anchor maps, and
clustering metadata when k-means anchors are used.

## Details

This is the seamless front-end advocated in the pair-observation
model-space RSA design: within-unit RSA fits and across-unit
representational connectivity become two views of the same fitted
object.

For `searchlight_result` inputs the full \\F F^\top\\ matrix is
intentionally never materialized. Instead, `k` anchor searchlights are
chosen (default: k-means clustering of the per-center fingerprint
matrix, with the searchlight closest to each centroid as the anchor) and
the connectivity is summarized as an `n_centers x k` similarity matrix
plus one brain map per anchor. Memory is \\\mathcal{O}(n\_{centers}
\cdot k)\\ rather than \\\mathcal{O}(n\_{centers}^2)\\.

## See also

[`rdm_model_space_connectivity`](http://bbuchsbaum.github.io/rMVPA/reference/rdm_model_space_connectivity.md),
[`pair_rsa_design`](http://bbuchsbaum.github.io/rMVPA/reference/pair_rsa_design.md),
[`rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md),
[`feature_rsa_rdm_vectors`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_rdm_vectors.md)
