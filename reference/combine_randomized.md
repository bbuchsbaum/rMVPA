# Combine randomized classifier results

This function combines the randomized classifier results from a good
results data frame and normalizes the performance matrix by the number
of instances for each voxel index.

## Usage

``` r
combine_randomized(model_spec, good_results, bad_results = NULL, ...)
```

## Arguments

- model_spec:

  A list containing the model specification.

- good_results:

  A data frame containing the successful classifier results.

- bad_results:

  A data frame containing the unsuccessful classifier results.

## Value

A list containing the combined and normalized performance matrix along
with other information from the dataset.
