# Create an MVPA Dataset Object

Creates a dataset object for MVPA analysis that encapsulates a training
dataset, an optional test dataset, and a voxel mask.

## Usage

``` r
mvpa_dataset(train_data, test_data = NULL, mask)
```

## Arguments

- train_data:

  The training data set: a `NeuroVec` instance

- test_data:

  Optional test data set: a `NeuroVec` instance (default: NULL)

- mask:

  The set of voxels to include: a `NeuroVol` instance

## Value

An `mvpa_dataset` object (S3 class) containing:

- train_data:

  The training data as a `NeuroVec` instance

- test_data:

  The test data as a `NeuroVec` instance (if provided, otherwise NULL)

- mask:

  The binary mask defining valid voxels as a `NeuroVol` instance

- has_test_set:

  Logical flag indicating whether this dataset has a test set

## See also

[`mvpa_surface_dataset`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_surface_dataset.md)
for creating surface-based MVPA datasets

[`mvpa_design`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
for creating the corresponding design object

## Examples

``` r
# Use gen_sample_dataset helper to create a simple dataset
sample_data <- gen_sample_dataset(c(5, 5, 5), nobs = 100, blocks = 4)
dataset <- sample_data$dataset

# Access components
print(dim(dataset$train_data))
#> [1]   5   5   5 100
print(sum(dataset$mask > 0))
#> [1] 125
```
