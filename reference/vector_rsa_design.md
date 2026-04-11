# Construct a design for a vectorized RSA model

This function constructs a design for an RSA model using a single
distance matrix, labels, and blocks.

## Usage

``` r
vector_rsa_design(D, labels, block_var)
```

## Arguments

- D:

  A representational dissimilarity matrix with row.names indicating
  labels.

- labels:

  character vector of labels corresponding to rows in another dataset X.

- block_var:

  A vector indicating the block (strata) each label belongs to. Must be
  the same length as \`labels\`.

## Value

A list containing the elements of the RSA design, class attributes
"vector_rsa_design" and "list".

## Details

The function verifies that all \`labels\` appear in \`rownames(D)\` and
creates an expanded dissimilarity matrix (\`Dexpanded\`) matching the
order of \`labels\`.

## Examples

``` r
# \donttest{
  D <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
  rownames(D) <- colnames(D) <- letters[1:5]
  labels <- factor(rep(letters[1:5], 2))
  block_var <- rep(1:2, each=5)
  des <- vector_rsa_design(D, labels, block_var)
# }
```
