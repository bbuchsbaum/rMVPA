# Evaluate model performance for feature RSA

Computes condition-pattern metrics (trial x trial correlation matrix),
voxel-level encoding metrics, global reconstruction metrics (MSE,
R-squared), and optionally performs permutation tests.

## Usage

``` r
evaluate_model.feature_rsa_model(
  object,
  predicted,
  observed,
  nperm = 0,
  save_distributions = FALSE,
  compute_rdm_vectors = isTRUE(object$return_rdm_vectors),
  fold_id = NULL,
  ...
)
```

## Arguments

- object:

  The feature RSA model

- predicted:

  Matrix of predicted values (observations x voxels)

- observed:

  Matrix of observed values (observations x voxels)

- nperm:

  Number of permutations for statistical testing (default: 0)

- save_distributions:

  Logical indicating whether to save full permutation distributions

- compute_rdm_vectors:

  Logical; when TRUE, also return compact predicted and observed RDM
  vectors for reuse by downstream code.

- fold_id:

  Optional vector assigning each observation to its outer
  cross-validation test fold. Pattern discrimination, identification
  rank, and RDM correlation then use only candidates or pairs withheld
  together. Cross-validated workflows supply this automatically.

- ...:

  Additional arguments

## Value

A list containing:

- pattern_correlation:

  Mean diagonal of the trial x trial correlation matrix – how well the
  predicted spatial pattern for each trial matches the correct observed
  pattern.

- pattern_discrimination:

  Diagonal minus off-diagonal of the trial x trial correlation matrix,
  restricted to candidates in the same outer test fold – how much better
  the correct trial is matched than eligible incorrect trials.

- pattern_rank_percentile:

  For each trial, percentile rank of the correct pattern among
  candidates in the same outer test fold. 0.5 = chance, 1 = perfect.

- rdm_correlation:

  Spearman correlation between predicted and observed
  correlation-distance pairs whose observations were withheld together.

- voxel_correlation:

  Correlation of the flattened predicted and observed matrices (global
  reconstruction quality).

- mse:

  Mean squared error.

- r_squared:

  1 - RSS/TSS.

- mean_voxelwise_temporal_cor:

  Average per-voxel temporal correlation (encoding fidelity).

- permutation_results:

  If `nperm > 0`, a list with p-values and z-scores for each metric.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Internal S3 method called after cross-validation
  # perf <- evaluate_model(feature_rsa_model, newdata, observed)
} # }
```
