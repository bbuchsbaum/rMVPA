# Perform resampled searchlight analysis

Similar to randomized searchlight but uses
neuroim2::resampled_searchlight, which draws a fixed total number of
centers (\`niter\`), optionally with a vector of radii. Results are
aggregated per voxel using the randomized combiner.

## Usage

``` r
do_resampled(
  model_spec,
  radius,
  niter,
  mvpa_fun = mvpa_iterate,
  combiner = combine_randomized,
  ...,
  drop_probs = FALSE,
  return_pobserved = FALSE,
  fail_fast = FALSE,
  backend = c("default", "shard", "auto")
)
```

## Arguments

- model_spec:

  An object specifying the model to be used in the searchlight analysis.

- radius:

  The radius (or vector of radii) of the searchlight sphere.

- niter:

  Total number of sampled searchlights.

- mvpa_fun:

  The MVPA function to be used in the searchlight analysis (default is
  `mvpa_iterate`).

- combiner:

  The function to be used to combine results (default is
  `combine_randomized`).

- ...:

  Additional arguments to be passed to the MVPA function.

- drop_probs:

  Logical; drop per-ROI probability matrices after computing metrics
  (default `FALSE`).

- return_pobserved:

  Logical; placeholder for API symmetry with randomized searchlight.
  Currently ignored because `combine_randomized` does not aggregate
  prob-observed.

## Value

A searchlight_result object containing spatial maps for each metric.
