# Fit ROI for contrast_rsa_model

Implements the `fit_roi` contract for MS-ReVE / contrast RSA models.
Called automatically by `process_roi.default` when a `fit_roi` method is
detected.

## Usage

``` r
# S3 method for class 'contrast_rsa_model'
fit_roi(model, roi_data, context, ...)
```

## Arguments

- model:

  A `contrast_rsa_model` object.

- roi_data:

  Named list with `train_data`, `indices`, `train_roi`, `test_roi`.

- context:

  Named list with `design`, `cv_spec`, `id`, `center_global_id`.

- ...:

  Additional arguments forwarded to `train_model`.

## Value

A `roi_result` object.
