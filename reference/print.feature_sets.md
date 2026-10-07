# Print method for feature_sets

Print method for feature_sets

## Usage

``` r
# S3 method for class 'feature_sets'
print(x, ...)
```

## Arguments

- x:

  feature_sets object

- ...:

  ignored

## Value

Invisibly returns the input object `x` (called for side effects).

## Examples

``` r
X <- matrix(rnorm(10 * 5), 10, 5)
fs <- feature_sets(X, blocks(low = 2, high = 3))
print(fs)
#> feature_sets
#> ===========
#> 
#> Observations: 10
#> Features:     5
#> Sets:         2 (low, high)
```
