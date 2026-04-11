# Get Center IDs for Searchlight Iteration

Returns the set of center IDs over which a searchlight analysis
iterates. For volumetric datasets this is the set of nonzero mask
voxels; for clustered datasets it is the sequence of cluster indices.

## Usage

``` r
get_center_ids(dataset, ...)
```

## Arguments

- dataset:

  The dataset object.

- ...:

  Additional arguments.

## Value

Integer vector of center IDs.

## Examples

``` r
ds <- gen_sample_dataset(c(5,5,5), 20)
ids <- get_center_ids(ds$dataset)
length(ids)
#> [1] 125
```
