# Apply Transformations to an Existing Contrast Matrix

Applies centering, scaling, and/or orthogonalization to a pre-existing
numeric contrast matrix.

## Usage

``` r
transform_contrasts(
  C,
  centre = TRUE,
  scale = c("none", "sd", "l2"),
  orth = FALSE,
  keep_attr = TRUE
)
```

## Arguments

- C:

  A numeric matrix representing the contrasts (conditions x contrasts).
  Row names, if present, should correspond to condition labels and will
  be preserved. Column names, if present, will be used for the
  \`"source"\` attribute if \`orth = TRUE\` and \`keep_attr = TRUE\`.

- centre:

  Logical. If TRUE (default), columns of the matrix are mean-centered.

- scale:

  Character string specifying scaling method after centering (if
  \`orth=FALSE\`). Options: \`"none"\` (default), \`"sd"\` (divide by
  sample standard deviation), \`"l2"\` (divide by L2 norm / vector
  length to get unit vectors). This argument is \*ignored\* if \`orth =
  TRUE\`.

- orth:

  Logical. If FALSE (default), the matrix columns represent the
  specified contrasts directly (after centering/scaling). If TRUE, an
  orthonormal basis for the column space is computed via QR
  decomposition. Resulting columns will be orthogonal and have unit
  length (L2 norm = 1).

- keep_attr:

  Logical. If TRUE (default) and \`orth = TRUE\`, the original column
  names (before orthogonalization) are stored in \`attr(C_transformed,
  "source")\`, and names of linearly dependent columns removed during
  orthogonalization are stored in \`attr(C_transformed, "dropped")\`.

## Value

A transformed numeric matrix. If \`orth = TRUE\` and \`keep_attr =
TRUE\`, it includes \`"source"\` and potentially \`"dropped"\`
attributes.

## Details

This function is useful for post-processing a contrast matrix,
especially one that might have been created by combining outputs from
different sources (e.g., theory-driven contrasts and data-driven
contrasts) or by direct manual construction.

## Scaling

Applied \*after\* centering if \`orth=FALSE\`.

- \`"none"\`: No scaling.

- \`"sd"\`: \`scale(..., center=FALSE, scale=TRUE)\`. Uses sample
  standard deviation (N-1 denominator). Note that for columns with few
  unique values (e.g., a centered +/-1 contrast), the SD can be slightly
  different depending on whether the number of items is even or odd, due
  to the N-1 denominator. This might lead to minor differences in scaled
  norms.

- \`"l2"\`: Divides each column by its L2 norm (\`sqrt(sum(x^2))\`).

## Orthogonalization

If \`orth = TRUE\`, uses \`qr.Q(qr(C))\` to find an orthonormal basis.
The number of columns in the output will be the rank of the input
matrix. Columns are renamed \`Orth1\`, \`Orth2\`, etc. Scaling is
ignored as the columns already have unit L2 norm. If \`keep_attr =
TRUE\`: \`attr(C_orth, "source")\` stores the names of the original
columns that formed the basis for the orthogonalized matrix.
\`attr(C_orth, "dropped")\` stores the names of original columns that
were linearly dependent and thus not part of the basis, if any.

## Specific Behaviors

- If \`orth = TRUE\` and the input matrix has only one column after
  potential centering, that column is scaled to unit L2 norm. Centering
  still depends on the \`centre\` argument.

- If \`centre = FALSE\` and \`orth = TRUE\`, the QR decomposition is
  performed on the \*uncentered\* columns.

- If the mini-DSL \`. \` notation is used for \`levelsB\` and
  \`levelsA\` already contains all \`labels\`, \`levelsB\` becomes
  empty, potentially resulting in a constant (zero) column before
  centering. A warning is issued in this case.

## See also

\[contrasts()\], \[make_feature_contrasts()\]

## Examples

``` r
C_manual <- matrix(c( 1, -1,
                      1,  1,
                      0,  0,
                      0,  0), nrow = 4, byrow = TRUE,
                   dimnames = list(paste0("Cond", 1:4), c("MainA", "MainB")))

# Center and make orthonormal
C_transformed <- transform_contrasts(C_manual, orth = TRUE)
print(C_transformed)
#>       Orth1         Orth2
#> Cond1  -0.5  7.071068e-01
#> Cond2  -0.5 -7.071068e-01
#> Cond3   0.5 -2.775558e-17
#> Cond4   0.5 -2.775558e-17
#> attr(,"source")
#> [1] "MainA" "MainB"
print(attr(C_transformed, "source"))
#> [1] "MainA" "MainB"

# Center and scale to unit L2 norm
C_l2 <- transform_contrasts(C_manual, scale = "l2")
print(round(colSums(C_l2^2), 5))
#> MainA MainB 
#>     1     1 
```
