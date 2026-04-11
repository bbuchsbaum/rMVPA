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
  [`orthogonalize_contrasts`](http://bbuchsbaum.github.io/rMVPA/reference/orthogonalize_contrasts.md).

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
if (FALSE) { # \dontrun{
# Example with non-overlapping contrasts (zero interaction)
C1 <- matrix(c(1,-1,0,0, 0,0,1,-1), nrow=4, 
             dimnames=list(NULL, c("A","B")))
# A compares conditions 1 vs 2, B compares 3 vs 4
# Their interaction will be zero and skipped

# Example with overlapping contrasts (non-zero interaction)  
C2 <- matrix(c(1,1,-1,-1, 1,-1,1,-1), nrow=4,
             dimnames=list(NULL, c("Main1","Main2")))
# These contrasts overlap and will produce a meaningful interaction
} # }
```
