# Compute a Neural RDM for RSA

Compute the full representational dissimilarity matrix (RDM) among
observed neural patterns using the correlation distances used by
[`rsa_model`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md).

## Usage

``` r
rsa_neural_rdm(
  patterns,
  method = c("spearman", "pearson"),
  pattern_center = c("none", "stimulus_mean")
)
```

## Arguments

- patterns:

  A finite numeric matrix with observations in rows and measured
  features (for example, voxels) in columns. At least two observations
  and two features are required. Each row must vary across features
  after optional pattern centering.

- method:

  Correlation method: `"spearman"` (the default, matching
  [`rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md))
  or `"pearson"`. Spearman uses average ranks for ties.

- pattern_center:

  Use `"none"` for uncentered patterns, or `"stimulus_mean"` to subtract
  each feature's mean across the supplied observations before computing
  correlations (and before Spearman ranking).

## Value

A symmetric numeric matrix with one row and column per observation, a
zero diagonal, and distances `1 - correlation`. Observation order and
row names are preserved on both axes; unnamed inputs yield unnamed
outputs.

## Details

These are correlation distances, not noise-whitened or crossvalidated
distances. No averaging, feature selection, run/item exclusion,
reliability estimate, or noise-ceiling calculation is performed. All
pairs are returned; exclusions belong to
[`rsa_design`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md)
or the caller.

To reproduce an
[`rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md)
neural RDM, supply the same preprocessed pattern matrix and feature
selection, with `method` matching the model's `distmethod`.
Stimulus-mean centering depends on all supplied observations: centering
separate runs or subsets generally changes the distances. The model
centers its input before any design-level row selection or pair
exclusions; reproduce that order when extracting subsets of this matrix.

Nonfinite inputs and patterns that are constant across features after
centering are rejected because their correlations are undefined.

## See also

[`rsa_model`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md),
[`rsa_design`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md),
[`cordist`](https://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)

## Examples

``` r
patterns <- rbind(item_a = c(1, 2, 4, 3),
                  item_b = c(2, 1, 3, 4),
                  item_c = c(4, 3, 1, 2))
rdm <- rsa_neural_rdm(patterns)
rdm[lower.tri(rdm)]
#> [1] 0.4 2.0 1.6
rsa_neural_rdm(patterns, method = "pearson",
               pattern_center = "stimulus_mean")
#>           item_a    item_b   item_c
#> item_a 0.0000000 0.6837722 1.857493
#> item_b 0.6837722 0.0000000 1.759257
#> item_c 1.8574929 1.7592566 0.000000
```
