# Perform randomized searchlight analysis

This function performs randomized searchlight analysis using a specified
model, radius, and number of iterations. It can be customized with
different MVPA functions, combiners, and permutation options.

## Usage

``` r
do_randomized(
  model_spec,
  radius,
  niter,
  mvpa_fun = mvpa_iterate,
  combiner = pool_randomized,
  ...,
  chunk_size = NULL,
  return_pobserved = TRUE,
  drop_probs = FALSE,
  fail_fast = FALSE,
  backend = c("default", "shard", "auto")
)
```

## Arguments

- model_spec:

  An object specifying the model to be used in the searchlight analysis.

- radius:

  The radius of the searchlight sphere.

- niter:

  The number of iterations for randomized searchlight.

- mvpa_fun:

  The MVPA function to be used in the searchlight analysis (default is
  `mvpa_iterate`).

- combiner:

  The function to be used to combine results (default is
  `pool_randomized`).

- ...:

  Additional arguments to be passed to the MVPA function.

## Value

A searchlight_result object containing spatial maps for each metric.
