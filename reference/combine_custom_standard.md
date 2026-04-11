# Combine Custom Standard Searchlight Results

Internal function to combine results from a standard custom searchlight
run. Creates a \`searchlight_result\` object with NeuroVol/NeuroSurface
for each metric.

## Usage

``` r
combine_custom_standard(dataset, iteration_results)
```

## Arguments

- dataset:

  The original mvpa_dataset object.

- iteration_results:

  The raw tibble output from \`mvpa_iterate\`.

## Value

A \`searchlight_result\` object.
