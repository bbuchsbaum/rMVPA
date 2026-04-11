# Extract Full Feature Matrix from a Dataset

Returns the full T x P feature matrix (observations x features) for
whole-brain analyses. For clustered datasets this is the parcel
time-series; for image datasets this is the voxel series under the mask.

## Usage

``` r
get_feature_matrix(dataset, ...)

# S3 method for class 'mvpa_clustered_dataset'
get_feature_matrix(dataset, ...)

# S3 method for class 'mvpa_image_dataset'
get_feature_matrix(dataset, ...)

# S3 method for class 'mvpa_multibasis_image_dataset'
get_feature_matrix(dataset, ...)

# S3 method for class 'mvpa_surface_dataset'
get_feature_matrix(dataset, ...)

# S3 method for class 'matrix'
get_feature_matrix(dataset, ...)
```

## Arguments

- dataset:

  An `mvpa_dataset` or a plain matrix.

- ...:

  Additional arguments passed to methods.

## Value

A numeric matrix of dimension T x P.

## Examples

``` r
ds <- gen_sample_dataset(c(5,5,5), 20)
X <- get_feature_matrix(ds$dataset)
dim(X)
#> [1]  20 125
```
