# Print Method for vector_rsa_design

Print Method for vector_rsa_design

## Usage

``` r
# S3 method for class 'vector_rsa_design'
print(x, ...)
```

## Arguments

- x:

  A vector_rsa_design object.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns the input object `x` (called for side effects).

## Examples

``` r
# \donttest{
  D <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
  rownames(D) <- colnames(D) <- letters[1:5]
  labels <- factor(rep(letters[1:5], 2))
  block_var <- rep(1:2, each=5)
  des <- vector_rsa_design(D, labels, block_var)
  print(des)
#> ================================================== \n         Vectorized RSA Design          \n================================================== \n\nDistance Matrix:\n  |- Original Dimensions:  5 x 5 \nLabels:\n  |- Total Labels:  10 \n  |- Sample:          a, b, c, d, e, ... \nBlocks:\n  |- Number of Blocks:  2 \n  |- Block Labels:      1, 2 \nExpanded D Matrix:\n  |- Dimensions:        10 x 10 \n\n ================================================== \n
# }
```
