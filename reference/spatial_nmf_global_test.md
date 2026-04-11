# Global Cross-validated Group Test for Spatial NMF

Tests whether groups are distinguishable in component space, without
attributing the effect to any single component. This is a multivariate,
cross-validated test: for each fold, NMF is fit on the training data and
component loadings for held-out subjects are used to train/test a
classifier. Permuting group labels yields a null distribution for the
chosen metric (AUC or accuracy).

## Usage

``` r
spatial_nmf_global_test(
  x = NULL,
  X = NULL,
  groups = NULL,
  k = NULL,
  lambda = NULL,
  graph = NULL,
  neighbors = 6,
  nfolds = 5,
  folds = NULL,
  nperm = 1000,
  permute = c("labels", "full"),
  metric = c("auc", "accuracy"),
  classifier = c("glm", "lda", "centroid"),
  scale = TRUE,
  positive = NULL,
  seed = NULL,
  project_args = list(max_iter = 100, check_every = 5),
  return_perm = FALSE,
  return_cv = FALSE,
  parallel = FALSE,
  future_seed = TRUE,
  progress = FALSE,
  ...
)
```

## Arguments

- x:

  Optional spatial_nmf_maps_result with return_data=TRUE.

- X:

  Optional subject-by-voxel matrix (overrides x\$data).

- groups:

  Factor or vector of group labels (length n).

- k:

  Number of components (defaults to x\$fit\$k if available).

- lambda:

  Spatial regularization strength (defaults to x\$fit\$lambda if
  available).

- graph:

  Optional graph Laplacian list (required if lambda \> 0). If
  \`graph\$A\` is weighted, set \`graph\$weighted=TRUE\` to preserve
  weights; otherwise edges are binarized.

- neighbors:

  Neighborhood size for volumetric adjacency (6/18/26).

- nfolds:

  Number of cross-validation folds.

- folds:

  Optional fold specification (vector of fold IDs or list of test
  indices). If NULL, folds are stratified by \`groups\`.

- nperm:

  Number of label permutations.

- permute:

  Permutation strategy: "labels" permutes labels with fixed folds and
  fixed NMF projections; "full" re-derives folds based on permuted
  labels and refits NMF within each fold. If \`folds\` is supplied,
  "full" falls back to "labels" (fixed folds).

- metric:

  Performance metric: "auc" or "accuracy".

- classifier:

  Classifier: "glm", "lda", or "centroid".

- scale:

  Logical; z-score W within each fold.

- positive:

  Optional positive class label (defaults to second factor level).

- seed:

  Optional RNG seed.

- project_args:

  List of arguments passed to spatial_nmf_project.

- return_perm:

  Logical; return permutation statistics.

- return_cv:

  Logical; return cross-validated predictions and fold IDs.

- parallel:

  Logical; use future_lapply for permutations (requires future.apply).

- future_seed:

  Optional seed control for future.apply (passed to future_lapply).

- progress:

  Logical; report progress via progressr (works with parallel futures).

- ...:

  Additional arguments passed to spatial_nmf_fit.

## Value

A list with the observed statistic, permutation p-value, and metadata.

## Details

Use this when you want a single omnibus test of group separability in
the learned component space.

- `stat` is the cross-validated performance (AUC or accuracy).

- `p_value` is the permutation p-value under the chosen strategy.

- `permute = "labels"` keeps folds fixed; `"full"` re-derives folds and
  refits NMF for each permutation (more expensive).

## Examples

``` r
if (FALSE) { # \dontrun{
  result <- spatial_nmf_global_test(
    matrix(rnorm(100*10), 100, 10),
    k = 3, nperm = 99
  )
} # }
```
