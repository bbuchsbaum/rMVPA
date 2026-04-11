# fit_roi method for era_rsa_model

Computes first-order ERA metrics and second-order encoding-retrieval
geometry for a single ROI/searchlight sphere.

## Usage

``` r
# S3 method for class 'era_rsa_model'
fit_roi(model, roi_data, context, ...)
```

## Arguments

- model:

  An era_rsa_model object.

- roi_data:

  List with train_data, test_data, indices, train_roi, test_roi.

- context:

  List with design, cv_spec, id, center_global_id.

- ...:

  Additional arguments (unused).

## Value

An `roi_result` object.
