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
#>              [,1]        [,2]
#>  [1,]  0.06640464  0.11308119
#>  [2,]  0.03448824 -0.02899078
#>  [3,]  0.05492971  0.29190164
#>  [4,]  0.17396037 -0.13989468
#>  [5,] -0.13015008 -0.03640639
#>  [6,]  0.01846279  0.10230708
#>  [7,]  0.13456863 -0.12180256
#>  [8,]  0.22313751 -0.09159190
#>  [9,] -0.07746549 -0.01223071
#> [10,] -0.06476447 -0.10092418
#> 
#> $importance
#>  [1] 0.13113707 0.04505445 0.29702498 0.22323247 0.13514610 0.10395967
#>  [7] 0.18150642 0.24120412 0.07842507 0.11991717
#> 
# }
```
