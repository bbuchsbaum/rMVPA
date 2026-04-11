# Declare the Output Metric Schema for a Model

Returns a named list describing the metrics that
[`fit_roi`](http://bbuchsbaum.github.io/rMVPA/reference/fit_roi.md)
produces for this model type. Used by the generic combiner (Phase 2) to
build the correct number of output maps without model-specific
`combine_*` functions.

## Usage

``` r
output_schema(model)

# S3 method for class 'manova_model'
output_schema(model)

# S3 method for class 'mvpa_model'
output_schema(model)

# S3 method for class 'item_model'
output_schema(model)

# S3 method for class 'rsa_model'
output_schema(model)

# S3 method for class 'vector_rsa_model'
output_schema(model)

# S3 method for class 'feature_rsa_model'
output_schema(model)

# S3 method for class 'hrfdecoder_model'
output_schema(model)

# S3 method for class 'naive_xdec_model'
output_schema(model)

# S3 method for class 'repmap_model'
output_schema(model)

# S3 method for class 'repmed_model'
output_schema(model)

# S3 method for class 'repnet_model'
output_schema(model)

# S3 method for class 'subspace_alignment_model'
output_schema(model)
```

## Arguments

- model:

  A model specification object.

## Value

A named list where names are metric names and values describe the type:
`"scalar"` (one value per ROI) or `"vector[N]"` (N values per ROI).
Returns `NULL` for models using the legacy combiner path.

## Examples

``` r
# \donttest{
  # Default returns NULL (legacy path)
  ds <- gen_sample_dataset(c(4,4,4), 20)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
  output_schema(mspec)  # NULL
#> $Accuracy
#> [1] "scalar"
#> 
#> $AUC
#> [1] "scalar"
#> 
# }
```
