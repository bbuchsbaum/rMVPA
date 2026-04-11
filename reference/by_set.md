# Define feature sets by a per-column set label

Use this when columns are not arranged as contiguous blocks.

## Usage

``` r
by_set(set, order = NULL)
```

## Arguments

- set:

  Character or factor vector of length ncol(X), naming the set for each
  column.

- order:

  Optional character vector giving the desired set order.

## Value

An object of class \`feature_set_spec_by\`.

## Details

This is useful when:

- you have interleaved predictors (e.g. time-lagged features),

- you want to group columns by an external annotation, or

- you already have column names encoding the set membership.

## Examples

``` r
X <- matrix(rnorm(10 * 6), 10, 6)
set <- rep(c("audio", "vision"), each = 3)
fs <- feature_sets(X, by_set(set, order = c("vision", "audio")))
fs
#> feature_sets
#> ===========
#> 
#> Observations: 10
#> Features:     6
#> Sets:         2 (vision, audio)
```
