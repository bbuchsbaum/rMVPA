# Number of Response Categories

Get the number of response categories or levels.

## Usage

``` r
nresponses(x)
```

## Arguments

- x:

  The object from which to extract the number of categories.

## Value

The number of response categories.

## Examples

``` r
des <- mvpa_design(
  train_design = data.frame(y = factor(c("a","b","a","b"))),
  y_train = factor(c("a","b","a","b")),
  block_var = factor(c("1","1","2","2"))
)
nresponses(des)
#> [1] 2
```
