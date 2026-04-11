# Format per-fold results: predict TR-level soft labels, aggregate to events

This method performs TR-level prediction on the test fold, then
aggregates the continuous predictions to event-level using the time
window specified in the model. Ground truth labels are extracted from
the events table during aggregation.

## Usage

``` r
# S3 method for class 'hrfdecoder_model'
format_result(obj, result, error_message = NULL, context, ...)
```

## Arguments

- obj:

  hrfdecoder_model specification object

- result:

  Output from train_model.hrfdecoder_model containing the fitted model

- error_message:

  Optional error message if training failed

- context:

  List containing test data and other fold information from crossval

- ...:

  Additional arguments (unused)

## Value

A tibble with one row containing class predictions, probabilities, true
labels, test indices, optional fit object, error status, and error
message

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires hrfdecoder package
  # See vignette("Continuous_Decoding") for full workflow
} # }
```
