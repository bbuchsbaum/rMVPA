# Orthogonalize a Contrast Matrix

Orthogonalizes the columns of a contrast matrix using QR decomposition.
The resulting matrix will have orthonormal columns spanning the same
space as the original columns, up to the rank of the input matrix. Sign
of the output columns is heuristically aligned with input columns.

## Usage

``` r
orthogonalize_contrasts(C)
```

## Arguments

- C:

  A numeric matrix (K x Q) where columns represent contrast vectors.

## Value

An orthogonalized matrix. If the input matrix `C` is rank-deficient
(rank \< Q), the output matrix will have Q columns, but only the first
`rank(C)` columns will be non-zero and form an orthonormal basis;
subsequent columns will be zero vectors. Column names from `C` are
preserved.

## Examples

``` r
K <- 6 # Number of conditions
Q <- 2 # Number of contrasts
C_orig <- matrix(c( 1,  1,  1, -1, -1, -1,  # Contrast 1
                    1, -1,  0,  1, -1,  0), # Contrast 2 (not orthogonal to C1)
                 nrow=K, ncol=Q)
colnames(C_orig) <- c("MainEffect", "InteractionLike")
C_ortho <- orthogonalize_contrasts(C_orig)
# print(round(crossprod(C_ortho), 10)) # Should be close to identity matrix

# Example with a rank-deficient matrix (3rd contrast is sum of first two)
C_rank_def <- cbind(C_orig, C_orig[,1] + C_orig[,2])
colnames(C_rank_def) <- c("C1", "C2", "C3_dependent")
C_ortho_def <- orthogonalize_contrasts(C_rank_def)
#> Warning: Input matrix C is rank-deficient (rank 2 < 3 columns). The output matrix has 3 columns, but only the first 2 are non-zero and form an orthonormal basis. Subsequent columns are zero vectors.
# print(round(crossprod(C_ortho_def), 10))
# The 3rd column of C_ortho_def will be zeros.
```
