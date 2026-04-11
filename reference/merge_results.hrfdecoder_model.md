# Merge fold results within an ROI to a classification_result and compute performance

This method concatenates event-level predictions across all CV folds,
filters out zero-probability events (which occur when events fall
outside the test fold), and builds a standard classification_result for
metric computation.

## Usage

``` r
# S3 method for class 'hrfdecoder_model'
merge_results(obj, result_set, indices, id, ...)
```

## Arguments

- obj:

  hrfdecoder_model specification object

- result_set:

  Tibble containing format_result outputs from all folds

- indices:

  ROI/searchlight voxel indices

- id:

  Unique identifier for this ROI/searchlight sphere

- ...:

  Additional arguments (unused)

## Value

A tibble with one row containing the classification_result, performance
metrics, optional fast metrics, primary metric name/value, and error
status

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires hrfdecoder package
  # See vignette("Continuous_Decoding") for full workflow
} # }
```
