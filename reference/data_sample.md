# Extract Sample from Dataset

Extract a sample from a given dataset object.

Construct a light-weight data sample descriptor for a set of features
(e.g., voxels) associated with a dataset. Methods define how to
interpret and convert this descriptor.

## Usage

``` r
data_sample(obj, vox, ...)

# S3 method for class 'mvpa_dataset'
data_sample(obj, vox, ...)

# S3 method for class 'mvpa_multibasis_image_dataset'
data_sample(obj, vox, ...)
```

## Arguments

- obj:

  An object (typically a dataset) from which to draw a sample.

- vox:

  The feature indices or coordinates defining the sample.

- ...:

  Additional arguments passed to methods.

## Value

A sample extracted from the dataset.

An object of class \`data_sample\` describing the sample; methods
provide conversions to ROI structures or data frames as needed.

## Examples

``` r
if (FALSE) { # \dontrun{
  ds <- gen_sample_dataset(c(5,5,5), 20)
  vox <- sample(which(ds$dataset$mask > 0), 10)
  samp <- data_sample(ds$dataset, vox)
} # }
```
