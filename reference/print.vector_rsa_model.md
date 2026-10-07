# Print Method for vector_rsa_model

Print Method for vector_rsa_model

## Usage

``` r
# S3 method for class 'vector_rsa_model'
print(x, ...)
```

## Arguments

- x:

  An object of class `vector_rsa_model`.

- ...:

  Additional arguments (ignored).

## Value

Invisibly returns the input object `x` (called for side effects).

## Examples

``` r
ds <- gen_sample_dataset(c(4, 4, 4), nobs = 10, blocks = 2)
D <- as.matrix(dist(matrix(rnorm(5 * 3), 5, 3)))
rownames(D) <- colnames(D) <- letters[1:5]
des <- vector_rsa_design(D, factor(rep(letters[1:5], 2)), rep(1:2, each = 5))
mdl <- vector_rsa_model(ds$dataset, des)
print(mdl)
#> = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#>           Vectorized RSA Model          
#> = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
#> 
#> Dataset:
#>   |- Data Dimensions:  4 x 4 x 4 x 10 
#>   |- Mask Length:      64 
#> 
#> Design:
#>   |- Number of Labels:     10 
#>   |- Number of Blocks:     2 
#>   |- Dissimilarity Matrix: 10 x 10 
#> 
#> Model Specification:
#>   |- Distance Function:    cordist 
#>   |- RSA Similarity Func:  pearson 
#> 
#> = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
```
