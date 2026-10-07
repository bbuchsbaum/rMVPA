# Combine MS-ReVE (Contrast RSA) Searchlight Results

This function gathers the Q-dimensional performance vectors from each
successful searchlight center and combines them into Q separate output
maps.

## Usage

``` r
combine_msreve_standard(model_spec, good_results, bad_results)
```

## Arguments

- model_spec:

  The `contrast_rsa_model` specification.

- good_results:

  A tibble containing successful results from
  `train_model.contrast_rsa_model`. Each row corresponds to a
  searchlight center. Expected columns include `id` (center voxel global
  index) and `performance` (a named numeric vector of length Q).

- bad_results:

  A tibble containing information about failed searchlights (for error
  reporting).

## Value

A `searchlight_result` object containing:

- results:

  A named list of `SparseNeuroVec` or `NeuroSurfaceVector` objects, one
  for each contrast (Q maps in total).

- ...:

  Other standard searchlight metadata.

## Examples

``` r
# Deprecated internal combiner; shown with mock per-center results
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 20, nlevels = 2, blocks = 2)
ids <- which(as.logical(ds$dataset$mask))[1:10]
good <- tibble::tibble(
  id = ids,
  performance = lapply(seq_along(ids), function(i) c(AvsB = rnorm(1), CvsD = rnorm(1)))
)
spec <- list(output_metric = "beta_delta", dataset = ds$dataset)
res <- suppressWarnings(combine_msreve_standard(spec, good, tibble::tibble()))
names(res$results)
#> [1] "AvsB" "CvsD"
```
