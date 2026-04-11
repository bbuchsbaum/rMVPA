# Combine Custom Randomized Searchlight Results

Internal function to combine results from a randomized custom
searchlight run. Averages results for each metric across overlapping
spheres.

## Usage

``` r
combine_custom_randomized(dataset, iteration_results)
```

## Arguments

- dataset:

  The original mvpa_dataset object.

- iteration_results:

  The raw tibble output from \*all\* iterations of \`mvpa_iterate\`.

## Value

A \`searchlight_result\` object.
