# Per-ROI REMAP fit

Implements the end-to-end REMAP flow for one ROI/searchlight: build
paired prototypes, whiten, fit reduced-rank mapping, transform source
trials, train base classifier in the target domain, and classify
target-domain trials. Emits a standard classification_result plus
ROI-level performance.

## Usage

``` r
# S3 method for class 'remap_rrr_model'
fit_roi(model, roi_data, context, ...)
```

## Value

An `roi_result` object.
