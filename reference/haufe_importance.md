# Haufe Feature Importance (Activation Patterns)

Computes per-feature importance using the Haufe et al. (2014)
transformation from decoding weights to encoding (activation) patterns:
`A = Sigma_x %*% W %*% solve(t(W) %*% Sigma_x %*% W)`.

## Usage

``` r
haufe_importance(
  W,
  Sigma_x = NULL,
  summary_fun = function(A) sqrt(rowSums(A^2)),
  X = NULL,
  center = TRUE,
  block_size = NULL
)
```

## Arguments

- W:

  A P x D weight matrix (features x discriminant directions).

- Sigma_x:

  A P x P covariance matrix of the training features. Ignored when `X`
  is supplied.

- summary_fun:

  A function applied to the rows of A to produce a scalar importance per
  feature. Defaults to the L2 norm across discriminants.

- X:

  Optional n x P matrix of training observations. When supplied, the
  matrix-free path is used and `Sigma_x` is not needed.

- center:

  Logical; centre the columns of `X` before computing the patterns
  (default `TRUE`, matching
  [`cov()`](https://rdrr.io/r/stats/cor.html)).

- block_size:

  Optional integer; when `X` is supplied, compute the activation
  patterns in column blocks of at most this many features.

## Value

A list with components:

- A:

  The P x D activation pattern matrix.

- importance:

  A numeric vector of length P with per-feature importance.

## Details

Two equivalent computation paths are available. Supplying `Sigma_x` uses
the explicit \\P \times P\\ covariance. Supplying the training
observations `X` instead uses the matrix-free identity \$\$Z = X_c W\$\$
\$\$A = X_c^\top Z (Z^\top Z)^{+}\$\$ where \\X_c\\ is the
column-centred data. The covariance normalisation cancels exactly, so
the two paths agree to numerical precision, but the matrix-free path
never forms a \\P \times P\\ matrix: it needs only \\n \times D\\, \\P
\times D\\, and \\D \times D\\ quantities. Use it for whole-brain
feature counts, where a dense covariance is infeasible. Because \\Z\\ is
centred, \\X_c^\top Z = X^\top Z\\, so `X` itself is never copied;
`block_size` additionally limits the number of feature columns
multiplied at once.

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
#>  [1,]  0.03640501 -0.175169969
#>  [2,]  0.07557471  0.013402844
#>  [3,] -0.11008962 -0.058702466
#>  [4,] -0.20928129  0.007006917
#>  [5,] -0.02694763 -0.002780939
#>  [6,] -0.03072565 -0.070472721
#>  [7,] -0.10072057  0.003351169
#>  [8,] -0.01765588 -0.077190825
#>  [9,] -0.11627070 -0.099868229
#> [10,]  0.08372862  0.147295206
#> 
#> $importance
#>  [1] 0.17891295 0.07675398 0.12476259 0.20939856 0.02709074 0.07687958
#>  [7] 0.10077631 0.07918430 0.15327276 0.16942951
#> 
  # identical, without forming cov(X):
  haufe_importance(W, X = X)
#> $A
#>              [,1]         [,2]
#>  [1,]  0.03640501 -0.175169969
#>  [2,]  0.07557471  0.013402844
#>  [3,] -0.11008962 -0.058702466
#>  [4,] -0.20928129  0.007006917
#>  [5,] -0.02694763 -0.002780939
#>  [6,] -0.03072565 -0.070472721
#>  [7,] -0.10072057  0.003351169
#>  [8,] -0.01765588 -0.077190825
#>  [9,] -0.11627070 -0.099868229
#> [10,]  0.08372862  0.147295206
#> 
#> $importance
#>  [1] 0.17891295 0.07675398 0.12476259 0.20939856 0.02709074 0.07687958
#>  [7] 0.10077631 0.07918430 0.15327276 0.16942951
#> 
# }
```
