# Get Unique Region IDs

Extract unique region IDs from a region mask, handling both volume and
surface data.

## Usage

``` r
get_unique_regions(region_mask, ...)
```

## Arguments

- region_mask:

  A region mask object (NeuroVol or NeuroSurface)

- ...:

  Additional arguments passed to methods

## Value

A sorted vector of unique region IDs greater than 0
