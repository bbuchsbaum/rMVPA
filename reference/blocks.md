# Define consecutive column blocks as feature sets

Convenience constructor for the common case where columns are arranged
as consecutive blocks (e.g., low/mid/high/sem PCs concatenated).

## Usage

``` r
blocks(...)
```

## Arguments

- ...:

  Named integers giving number of columns per set (must sum to
  \`ncol(X)\`).

## Value

An object of class \`feature_set_spec_blocks\`.

## Details

\`blocks()\` does not touch any data; it only declares how columns
should be interpreted when passed to \`feature_sets()\`.

## Examples

``` r
spec <- blocks(low = 100, mid = 100, high = 100, sem = 100)
```
