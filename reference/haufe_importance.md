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
#>              [,1]        [,2]
#>  [1,] -0.20601881 -0.04836532
#>  [2,] -0.26969162 -0.07218315
#>  [3,]  0.02941763 -0.15261698
#>  [4,] -0.06307709  0.04272807
#>  [5,]  0.14622583 -0.27530518
#>  [6,] -0.06839436 -0.01409065
#>  [7,] -0.14964591  0.13610741
#>  [8,]  0.09851969 -0.19713427
#>  [9,] -0.25944046 -0.16014209
#> [10,]  0.17747513  0.08697459
#> 
#> $importance
#>  [1] 0.21161984 0.27918448 0.15542632 0.07618666 0.31172895 0.06983075
#>  [7] 0.20228476 0.22038159 0.30488496 0.19764109
#> 
  # identical, without forming cov(X):
  haufe_importance(W, X = X)
#> $A
#>              [,1]        [,2]
#>  [1,] -0.20601881 -0.04836532
#>  [2,] -0.26969162 -0.07218315
#>  [3,]  0.02941763 -0.15261698
#>  [4,] -0.06307709  0.04272807
#>  [5,]  0.14622583 -0.27530518
#>  [6,] -0.06839436 -0.01409065
#>  [7,] -0.14964591  0.13610741
#>  [8,]  0.09851969 -0.19713427
#>  [9,] -0.25944046 -0.16014209
#> [10,]  0.17747513  0.08697459
#> 
#> $importance
#>  [1] 0.21161984 0.27918448 0.15542632 0.07618666 0.31172895 0.06983075
#>  [7] 0.20228476 0.22038159 0.30488496 0.19764109
#> 
# }
```
