# Print Method for Feature RSA Design

Print Method for Feature RSA Design

## Usage

``` r
# S3 method for class 'feature_rsa_design'
print(x, ...)
```

## Arguments

- x:

  A feature_rsa_design object.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns the input object `x` (called for side effects).

## Examples

``` r
# \donttest{
  S <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
  labels <- factor(letters[1:5])
  des <- feature_rsa_design(S = S, labels = labels)
  print(des)
#> ================================================== 
#>           Feature RSA Design          
#> ================================================== 
#> 
#> Number of Observations:  5 
#> Feature Dimensions:      2 
#> Max Components Limit:    2 
#> Blocking Variable:       None
#> Similarity Matrix:       Provided
#> Labels (first few):    a, b, c, d, e 
#> 
#>  ================================================== 
# }
```
