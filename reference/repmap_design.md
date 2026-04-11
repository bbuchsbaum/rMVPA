# Representational mapping design helper (ReNA-Map)

Construct a design object for representational mapping, specifying
item-wise seed feature vectors.

## Usage

``` r
repmap_design(items, seed_features)
```

## Arguments

- items:

  Character vector of item IDs (keys).

- seed_features:

  K x P matrix of seed features; rows = items.

## Value

A list with elements:

- `items`: character vector of items

- `seed_features`: K x P feature matrix with rownames = items

## Examples

``` r
if (FALSE) { # \dontrun{
  items <- letters[1:10]
  seed_features <- matrix(rnorm(10*5), 10, 5)
  des <- repmap_design(items, seed_features)
} # }
```
