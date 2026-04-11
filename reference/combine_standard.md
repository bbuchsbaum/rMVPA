# Combine standard classifier results

This function combines the standard classifier results from a good
results data frame by binding the performance rows together and
optionally computes the observed probabilities.

## Usage

``` r
combine_standard(model_spec, good_results, bad_results)
```

## Arguments

- model_spec:

  A list containing the model specification

- good_results:

  A data frame containing the successful classifier results

- bad_results:

  A data frame containing the unsuccessful classifier results

## Value

A list containing the combined performance matrix and other information
