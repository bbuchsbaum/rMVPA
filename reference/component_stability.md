# Compare spatial subspaces across retained folds

Compare spatial subspaces across retained folds

## Usage

``` r
component_stability(object)
```

## Arguments

- object:

  A global result retaining at least two fits.

## Value

One row per fold pair, with effective ranks, number of common retained
features, principal angles in radians (a list column), and
`overlap = sum(cos(angles)^2) / max(rank1, rank2)`. The normalization
penalizes unequal ranks. Empty subspaces have undefined (NA) overlap.

## Details

Comparison uses original feature units on the intersection of columns
retained by both folds. It is invariant to orthogonal component
rotations. It measures spatial subspaces, not matched components or
statistical evidence for a supported rank.
