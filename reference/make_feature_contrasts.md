# Generate Contrasts from a Feature Matrix (Optional PCA)

Creates contrasts based on a matrix where rows represent conditions and
columns represent features (e.g., neural network embeddings, semantic
features). Optionally performs PCA to reduce dimensionality.

## Usage

``` r
make_feature_contrasts(
  features,
  labels = NULL,
  use_pca = TRUE,
  centre_pca = TRUE,
  scale_pca = FALSE,
  pve = 0.9,
  n_pcs = NULL,
  prefix = "Feat_"
)
```

## Arguments

- features:

  A numeric matrix (K x P) where K is the number of conditions and P is
  the number of features. Row names, if present, should correspond to
  condition labels. Column names are recommended.

- labels:

  Optional character vector of condition labels. If provided, rows of
  \`features\` matrix will be reordered to match this order. If NULL,
  the order from \`rownames(features)\` is used (if available).

- use_pca:

  Logical. If TRUE (default), performs Principal Component Analysis
  (PCA) on the features. If FALSE, uses the raw features directly.

- centre_pca:

  Logical. If \`use_pca = TRUE\`, should features be centered before
  PCA? (Default: TRUE) Note: Reordering of rows based on \`labels\`
  argument happens \*before\* PCA.

- scale_pca:

  Logical. If \`use_pca = TRUE\`, should features be scaled to unit
  variance before PCA? (Default: FALSE, as scaling can affect variance
  explained). Note: Reordering of rows based on \`labels\` argument
  happens \*before\* PCA.

- pve:

  Numeric (0 to 1). If \`use_pca = TRUE\`, selects the minimum number of
  principal components (PCs) needed to explain at least this proportion
  of variance. Ignored if \`n_pcs\` is specified. (Default: 0.9)

- n_pcs:

  Integer. If \`use_pca = TRUE\`, selects exactly this number of
  principal components. Takes precedence over \`pve\`. (Default: NULL)

- prefix:

  Character string to prepend to column names of the output matrix
  (e.g., "Feat\_", "PCA\_"). (Default: "Feat\_")

## Value

A numeric matrix (K x Q) where K matches the number of conditions/labels
and Q is the number of selected features or principal components. Rows
are ordered according to \`labels\` or \`rownames(features)\`. Columns
are named using the \`prefix\` and either the original feature names (if
\`use_pca=FALSE\`) or component numbers (e.g., "PCA_PC1", "PCA_PC2").

## See also

\[contrasts()\], \[transform_contrasts()\]

## Examples

``` r
# Example feature matrix (4 conditions, 5 features)
feat_mat <- matrix(rnorm(20), nrow = 4,
                   dimnames = list(paste0("Cond", 1:4), paste0("F", 1:5)))

# Use raw features (first 3)
C_raw <- make_feature_contrasts(feat_mat[, 1:3], use_pca = FALSE, prefix="RawFeat_")
print(C_raw)
#>        RawFeat_F1  RawFeat_F2  RawFeat_F3
#> Cond1 -2.92611970 -0.01410577  0.40531860
#> Cond2  0.02866434  0.81715471  1.74367917
#> Cond3  1.58705327  1.27925049 -0.05634012
#> Cond4 -0.43410305  0.51550743 -0.03770414

# Use PCA, selecting top 2 PCs
C_pca <- make_feature_contrasts(feat_mat, use_pca = TRUE, n_pcs = 2, prefix="PCA_")
print(C_pca)
#>          PCA_PC1    PCA_PC2
#> Cond1  2.7794656 -0.8641491
#> Cond2 -0.6814743 -0.4630886
#> Cond3 -3.2361019 -0.1316229
#> Cond4  1.1381106  1.4588605

# Use PCA, selecting >= 80% variance explained
C_pca_pve <- make_feature_contrasts(feat_mat, use_pca = TRUE, pve = 0.8, prefix="PCA_")
print(C_pca_pve)
#>          PCA_PC1    PCA_PC2
#> Cond1  2.7794656 -0.8641491
#> Cond2 -0.6814743 -0.4630886
#> Cond3 -3.2361019 -0.1316229
#> Cond4  1.1381106  1.4588605

# Reorder based on labels
C_pca_reorder <- make_feature_contrasts(feat_mat, labels=c("Cond3", "Cond1", "Cond4", "Cond2"),
                                      use_pca = TRUE, n_pcs = 2, prefix="PCA_")
print(C_pca_reorder)
#>          PCA_PC1    PCA_PC2
#> Cond3  3.2361019  0.1316229
#> Cond1 -2.7794656  0.8641491
#> Cond4 -1.1381106 -1.4588605
#> Cond2  0.6814743  0.4630886
```
