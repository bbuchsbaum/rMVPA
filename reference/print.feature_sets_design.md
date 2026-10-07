# Print method for feature_sets_design

Print method for feature_sets_design

## Usage

``` r
# S3 method for class 'feature_sets_design'
print(x, ...)
```

## Arguments

- x:

  A feature_sets_design object

- ...:

  ignored

## Value

Invisibly returns the input object `x` (called for side effects).

## Examples

``` r
X <- matrix(rnorm(10 * 5), 10, 5)
fs <- feature_sets(X, blocks(low = 2, high = 3))
des <- feature_sets_design(fs, block_var_train = rep(1:2, each = 5))
print(des)
#> feature_sets_design
#> =================
#> 
#> Train (encoding): 10 x 5 (sets: low, high)
#> Test (recall):    <none>
```
