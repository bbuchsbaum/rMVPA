# Compute performance for hrfdecoder_model

This method calls the performance function stored in the model
specification to compute metrics on the classification_result produced
by merge_results.

## Usage

``` r
# S3 method for class 'hrfdecoder_model'
compute_performance(obj, result)
```

## Arguments

- obj:

  hrfdecoder_model specification object

- result:

  A classification_result object from merge_results

## Value

Named numeric vector of performance metrics (e.g., Accuracy, AUC)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires hrfdecoder package
  # See vignette("Continuous_Decoding") for full workflow
} # }
```
