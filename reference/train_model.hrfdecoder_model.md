# Train per ROI/fold using hrfdecoder

This method trains the continuous-time decoder on TR-level data from the
training fold. Unlike standard MVPA, it does not use the \`y\` parameter
for training - the actual decoding targets come from the event model and
events stored in the design object.

## Usage

``` r
# S3 method for class 'hrfdecoder_model'
train_model(obj, train_dat, y, sl_info, cv_spec, indices, ...)
```

## Arguments

- obj:

  hrfdecoder_model spec

- train_dat:

  Tibble/matrix of TR x features (ROI data for train rows)

- y:

  IGNORED. This is a dummy TR sequence from CV fold construction. Actual
  targets come from obj\$design\$event_model and obj\$design\$events.

- sl_info:

  Searchlight info list with center ids (optional)

- cv_spec:

  Cross-validation spec (optional)

- indices:

  Global ROI indices

- ...:

  Additional arguments (unused)

## Value

A list with class "hrfdecoder_fit_wrap" containing the fitted model,
searchlight info, and ROI indices

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires hrfdecoder package
  # See vignette("Continuous_Decoding") for full workflow
} # }
```
