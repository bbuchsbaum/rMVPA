# Merge Results for MANOVA Model

This function takes the computed -log(p-values) from
\`train_model.manova_model\` for a single ROI/searchlight and formats it
into the standard output tibble.

This function takes the computed coefficients/correlations/t-values from
`train_model.rsa_model` for a single ROI/searchlight and formats it into
the standard output tibble.

Aggregates results (scores) and calls evaluate_model. Vector RSA
typically doesn't involve folds in the same way as classifiers, so this
mainly formats the output of train_model for the specific ROI/sphere.

## Usage

``` r
# S3 method for class 'manova_model'
merge_results(obj, result_set, indices, id, ...)

# S3 method for class 'mvpa_model'
merge_results(obj, result_set, indices, id, ...)

# S3 method for class 'binary_classification_result'
merge_results(obj, ...)

# S3 method for class 'regression_result'
merge_results(obj, ...)

# S3 method for class 'multiway_classification_result'
merge_results(obj, ...)

# S3 method for class 'regional_mvpa_result'
merge_results(obj, ...)

# S3 method for class 'rsa_model'
merge_results(obj, result_set, indices, id, ...)

# S3 method for class 'vector_rsa_model'
merge_results(obj, result_set, indices, id, ...)

# S3 method for class 'feature_rsa_model'
merge_results(obj, result_set, indices, id, ...)
```

## Arguments

- obj:

  The vector RSA model specification (contains nperm etc.).

- result_set:

  A tibble from the processor. Expected to contain the output of
  \`train_model.vector_rsa_model\` (the scores vector) likely within
  \`\$result\[\[1\]\]\`.

- indices:

  Voxel indices for the current ROI/searchlight sphere.

- id:

  Identifier for the current ROI/searchlight center.

- ...:

  Additional arguments.

## Value

A tibble row with the formatted performance metrics for the ROI/sphere.

A merged `regional_mvpa_result` object.

A tibble row with the formatted "performance" metrics
(coefficients/t-values/correlations) for the ROI/sphere.

A tibble row with the final performance metrics for the ROI/sphere.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Internal S3 method called during ROI processing
  # result <- merge_results(manova_model, result_set, indices, id)
} # }
```
