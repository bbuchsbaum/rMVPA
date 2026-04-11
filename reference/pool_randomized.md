# Combine randomized searchlight results by pooling

This function combines randomized searchlight results by pooling the
good results.

## Usage

``` r
pool_randomized(
  model_spec,
  good_results,
  bad_results = NULL,
  chunk_size = NULL,
  return_pobserved = TRUE,
  ...
)
```

## Arguments

- model_spec:

  An object specifying the model used in the searchlight analysis.

- good_results:

  A data frame containing the valid searchlight results.

- bad_results:

  A data frame containing the invalid searchlight results.

## Value

An object containing the combined searchlight results.
