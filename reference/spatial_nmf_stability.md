# Bootstrap Stability for Spatial NMF Components

Quantifies how stable the learned component maps are to resampling
subjects. For each bootstrap (or subsample), the NMF is re-fit,
components are matched to the reference solution, and summary maps are
accumulated.

## Usage

``` r
spatial_nmf_stability(
  x = NULL,
  X = NULL,
  fit = NULL,
  graph = NULL,
  lambda = NULL,
  n_boot = 200,
  sample = c("bootstrap", "subsample"),
  sample_frac = 1,
  init = c("nndsvd", "random"),
  fast = FALSE,
  normalize = c("H", "none"),
  similarity = c("cosine", "cor"),
  match = c("greedy"),
  top_frac = 0.1,
  seed = NULL,
  return_maps = FALSE,
  parallel = FALSE,
  future_seed = TRUE,
  progress = FALSE,
  ...
)
```

## Arguments

- x:

  A spatial_nmf_maps_result or spatial_nmf_fit object.

- X:

  Optional data matrix (n x p) if not included in x.

- fit:

  Optional spatial_nmf_fit object (if x is not provided).

- graph:

  Optional graph Laplacian list (required if lambda \> 0).

- lambda:

  Spatial regularization strength (defaults to fit\$lambda).

- n_boot:

  Number of bootstrap samples.

- sample:

  One of "bootstrap" or "subsample".

- sample_frac:

  Fraction of subjects to sample.

- init:

  NMF initialization for bootstrap fits.

- fast:

  Logical; use faster defaults for each bootstrap fit (e.g., random
  init, fewer iterations). You can still override specific optimization
  settings via \`...\`.

- normalize:

  Component normalization ("H" rescales rows to sum 1).

- similarity:

  Similarity measure for component matching ("cosine" or "cor").

- match:

  Matching strategy (currently "greedy").

- top_frac:

  Fraction of top voxels used to compute selection frequency.

- seed:

  Optional RNG seed.

- return_maps:

  Logical; return stability maps as NeuroVol/NeuroSurface.

- parallel:

  Logical; use future_lapply for bootstrap resamples (requires
  future.apply).

- future_seed:

  Optional seed control for future.apply (passed to future_lapply).

- progress:

  Logical; report progress via progressr (works with parallel futures).

- ...:

  Additional arguments passed to spatial_nmf_fit.

## Value

A list with stability summaries and optional maps.

## Details

Use this to assess whether components are reproducible and which voxels
are consistently among the strongest loadings.

- `mean`, `sd`, `cv`: bootstrap mean/SD/CV of component maps.

- `selection`: frequency with which a voxel appears in the top fraction.

- `component_similarity`: average similarity to the reference
  components.

## Examples

``` r
if (FALSE) { # \dontrun{
  stab <- spatial_nmf_stability(
    matrix(rnorm(100*10), 100, 10),
    k = 3, nruns = 5
  )
} # }
```
