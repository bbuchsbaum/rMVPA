# Add Interaction Contrasts to an msreve_design

Creates new contrast columns representing pairwise interactions of
existing contrasts in an `msreve_design` object. Interactions are
computed as element-wise products of the contrast vectors.

## Usage

``` r
add_interaction_contrasts(design, pairs = NULL, orthogonalize = TRUE)
```

## Arguments

- design:

  An object of class `msreve_design`.

- pairs:

  Optional two-column matrix or list of character vectors specifying
  pairs of contrast column names. Default `NULL` uses all pairwise
  combinations.

- orthogonalize:

  Logical; if `TRUE` (default) the expanded contrast matrix is passed
  through
  [`orthogonalize_contrasts`](https://bbuchsbaum.github.io/rMVPA/reference/orthogonalize_contrasts.md).

## Value

The updated `msreve_design` object with non-zero interaction columns
appended. Zero interactions are automatically skipped.

## Details

Interaction contrasts are created by element-wise multiplication of
pairs of contrast vectors. If the resulting interaction is a zero vector
(which occurs when the original contrasts have non-overlapping support,
i.e., no conditions where both contrasts are non-zero), the interaction
is skipped with an informative message. This commonly happens with
contrasts that compare distinct subsets of conditions, such as
c(1,-1,0,0) and c(0,0,1,-1).

## Examples

``` r
# Design with 4 conditions
df <- data.frame(cond = factor(rep(c("c1", "c2", "c3", "c4"), times = 6)),
                 run  = rep(1:3, each = 8))
mvdes <- mvpa_design(df, y_train = ~ cond, block_var = ~ run)

# Overlapping contrasts: their interaction is non-zero and is added
C2 <- matrix(c(1, 1, -1, -1,
               1, -1, 1, -1), nrow = 4,
             dimnames = list(c("c1", "c2", "c3", "c4"), c("Main1", "Main2")))
des2 <- add_interaction_contrasts(msreve_design(mvdes, C2))
colnames(des2$contrast_matrix)
#> [1] "Main1"         "Main2"         "Main1_x_Main2"

# Non-overlapping contrasts (1 vs 2, 3 vs 4): interaction is zero and skipped
C1 <- matrix(c(1, -1, 0, 0,
               0, 0, 1, -1), nrow = 4,
             dimnames = list(c("c1", "c2", "c3", "c4"), c("A", "B")))
des1 <- add_interaction_contrasts(msreve_design(mvdes, C1))
#> Interaction A_x_B is zero (contrasts have non-overlapping support) and will be skipped
colnames(des1$contrast_matrix)
#> [1] "A" "B"
```
