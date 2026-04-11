# Create a Feature-Based RSA Design

Creates a design for feature-based Representational Similarity Analysis
(RSA). You can either supply a similarity matrix S (and optionally
select dimensions) or directly supply a feature matrix F.

## Usage

``` r
feature_rsa_design(
  S = NULL,
  F = NULL,
  labels,
  k = 0,
  max_comps = 10,
  block_var = NULL
)
```

## Arguments

- S:

  A symmetric similarity matrix representing the feature space
  relationships. If NULL, you must supply F.

- F:

  A feature space matrix (observations by features). If supplied, this
  overrides S and k.

- labels:

  Vector of labels corresponding to the rows/columns of S or
  observations of F.

- k:

  Integer specifying the number of feature dimensions to retain when
  using S. If 0 (default), automatically determines dimensions using
  eigenvalue threshold \> 1 (minimum 2 dimensions kept). This parameter
  is ignored if F is supplied directly (k becomes ncol(F)).

- max_comps:

  Initial upper limit for the number of components to be derived from
  the feature space F by subsequent \`feature_rsa_model\` methods (PCA,
  PLS). This value is automatically capped by the final feature
  dimensionality \`k\`. Default 10.

- block_var:

  Optional blocking variable for cross-validation. If provided and
  \`crossval\` is \`NULL\` in \`feature_rsa_model\`, a blocked
  cross-validation scheme will be generated using this vector.

## Value

A `feature_rsa_design` object (S3 class) containing:

- S:

  The input similarity matrix (if used)

- F:

  Feature space projection matrix (k dimensions)

- labels:

  Vector of observation labels

- k:

  The final number of feature dimensions used

- max_comps:

  The upper limit on components (\<= k)

- block_var:

  Optional blocking variable for cross-validation

## Details

This function defines the feature space representation for the analysis.
If F is supplied directly, it is used as-is, and \`k\` becomes
\`ncol(F)\`. If only S is supplied, an eigen decomposition of S is
performed. \`k\` determines how many eigenvectors form the feature
matrix F. If \`k=0\`, dimensions with eigenvalues \> 1 are kept (minimum
2). \`max_comps\` sets an upper bound for the number of components that
model-fitting methods (like PCA, PLS in \`feature_rsa_model\`) can use,
and it cannot exceed the final feature dimensionality \`k\`.

## Examples

``` r
# \donttest{
  S <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
  labels <- factor(letters[1:5])
  des <- feature_rsa_design(S = S, labels = labels)
# }
```
