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
if (FALSE) { # \dontrun{
  # Internal function for combining MS-ReVE results
  # result <- combine_msreve_standard(model_spec, good_results, bad_results)
} # }
```
