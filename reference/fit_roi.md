# Fit a Model on a Single ROI

The primary dispatch point for per-ROI analysis in the new architecture.
Each model type implements this method to encapsulate its complete
analysis pipeline (including any internal cross-validation).

Computes mediation paths a, b, c' and the indirect effect in RDM space
for a single ROI/searchlight.

Computes representational connectivity metrics for a single
ROI/searchlight.

## Usage

``` r
fit_roi(model, roi_data, context, ...)

# S3 method for class 'manova_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'mvpa_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'item_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'rsa_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'vector_rsa_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'feature_rsa_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'hrfdecoder_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'naive_xdec_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'repmap_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'repmed_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'repnet_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'subspace_alignment_model'
fit_roi(model, roi_data, context, ...)

# S3 method for class 'feature_rsa_da_model'
fit_roi(model, roi_data, context, ...)
```

## Arguments

- model:

  The model specification object.

- roi_data:

  A list with components:

  train_data

  :   Numeric matrix (observations x features)

  test_data

  :   Numeric matrix or NULL

  indices

  :   Integer vector of voxel/vertex indices

  train_roi

  :   The raw ROI object (for models needing spatial info)

  test_roi

  :   The raw test ROI object or NULL

- context:

  A list with components:

  design

  :   The design object

  cv_spec

  :   Cross-validation specification or NULL

  id

  :   ROI identifier (center voxel ID or region number)

  center_global_id

  :   Global index of center voxel or NA

- ...:

  Additional model-specific arguments.

## Value

A
[`roi_result`](http://bbuchsbaum.github.io/rMVPA/reference/roi_result.md)
object.

## Details

Models that implement `fit_roi` are automatically preferred by
[`process_roi.default`](http://bbuchsbaum.github.io/rMVPA/reference/process_roi-methods.md)
over the legacy dispatch chain (`internal_crossval` /
`external_crossval` / `train_model`).

## Architecture TODO

`fit_roi` is currently scoped to ROI/searchlight iteration. Whole-brain
global analysis uses
[`run_global`](http://bbuchsbaum.github.io/rMVPA/reference/run_global.md)
via a separate `cv_run_global -> train_model` path. A future refactor
may unify these flows under a single fit contract.

## See also

[`roi_result`](http://bbuchsbaum.github.io/rMVPA/reference/roi_result.md),
[`output_schema`](http://bbuchsbaum.github.io/rMVPA/reference/output_schema.md),
[`process_roi`](http://bbuchsbaum.github.io/rMVPA/reference/process_roi-methods.md)

## Examples

``` r
# \donttest{
  # fit_roi is typically called internally by process_roi.default.
  # To implement for a new model class:
  # fit_roi.my_model <- function(model, roi_data, context, ...) {
  #   metrics <- c(accuracy = 0.85)
  #   roi_result(metrics, roi_data$indices, context$id)
  # }
# }
```
