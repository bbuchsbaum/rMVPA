# Get Multiple Data Samples

Extract multiple data samples based on a list of voxel/index sets from a
dataset object.

## Usage

``` r
get_samples(obj, vox_list)
```

## Arguments

- obj:

  The input dataset object.

- vox_list:

  A list of vectors containing voxel indices to extract.

## Value

A list of data samples.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks=2)
  vox_list <- list(1:10, 11:20)
  samples <- get_samples(ds$dataset, vox_list)
# }
```
