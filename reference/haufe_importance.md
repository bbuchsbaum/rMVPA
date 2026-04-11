# Haufe Feature Importance (Activation Patterns)

Computes per-feature importance using the Haufe et al. (2014)
transformation from decoding weights to encoding (activation) patterns:
`A = Sigma_x %*% W %*% solve(t(W) %*% Sigma_x %*% W)`.

## Usage

``` r
haufe_importance(W, Sigma_x, summary_fun = function(A) sqrt(rowSums(A^2)))
```

## Arguments

- W:

  A P x D weight matrix (features x discriminant directions).

- Sigma_x:

  A P x P covariance matrix of the training features.

- summary_fun:

  A function applied to the rows of A to produce a scalar importance per
  feature. Defaults to the L2 norm across discriminants.

## Value

A list with components:

- A:

  The P x D activation pattern matrix.

- importance:

  A numeric vector of length P with per-feature importance.

## References

Haufe, S., Meinecke, F., Goergen, K., Daehne, S., Haynes, J.D.,
Blankertz, B., & Biessmann, F. (2014). On the interpretation of weight
vectors of linear models in multivariate neuroimaging. NeuroImage, 87,
96-110.

## Examples

``` r
# \donttest{
  W <- matrix(rnorm(10*2), 10, 2)
  X <- matrix(rnorm(50*10), 50, 10)
  haufe_importance(W, cov(X))
#> $A
#>              [,1]         [,2]
#>  [1,] -0.18550545 -0.088476345
#>  [2,] -0.16381826 -0.033305381
#>  [3,]  0.04595221 -0.020781199
#>  [4,] -0.10096699 -0.111831456
#>  [5,]  0.13366036 -0.038591991
#>  [6,] -0.11836237  0.158550640
#>  [7,] -0.06055844  0.238514958
#>  [8,]  0.03028341 -0.166105527
#>  [9,] -0.01724757 -0.017244576
#> [10,]  0.04042202  0.001290194
#> 
#> $importance
#>  [1] 0.20552453 0.16716959 0.05043277 0.15066721 0.13912022 0.19785843
#>  [7] 0.24608273 0.16884351 0.02438963 0.04044261
#> 
# }
```
