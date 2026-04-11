# Component-level Inference for Spatial NMF

Tests whether subjects express NMF components differently between
groups, using the component loadings \\W\\ from a spatial NMF fit. This
is a component-wise analysis (not voxelwise): each component gets a
statistic and a permutation p-value.

## Usage

``` r
spatial_nmf_component_test(
  fit = NULL,
  W = NULL,
  groups = NULL,
  covariates = NULL,
  test = c("two_group", "one_group"),
  nperm = 1000,
  correction = c("maxT", "none"),
  alternative = c("two.sided", "greater"),
  null_W = NULL,
  alpha = 0.05,
  seed = NULL,
  return_perm = FALSE,
  parallel = FALSE,
  future_seed = TRUE,
  progress = FALSE
)
```

## Arguments

- fit:

  A spatial_nmf_fit object containing a W matrix.

- W:

  Optional n x k matrix of loadings (overrides fit\$W).

- groups:

  Factor or vector of group labels (length n); required for two-group
  tests.

- covariates:

  Optional data frame of covariates (n rows).

- test:

  One of "two_group" or "one_group".

- nperm:

  Number of permutations (ignored for one_group if null_W provided).

- correction:

  One of "maxT" or "none".

- alternative:

  Alternative hypothesis: "two.sided" or "greater". For two_group tests,
  "greater" evaluates positive differences for the second factor level
  relative to the first.

- null_W:

  Null distribution of W for one-group inference (list or 3D array).

- alpha:

  Significance threshold for counting significant components.

- seed:

  Optional RNG seed for permutations.

- return_perm:

  Logical; return permutation statistics.

- parallel:

  Logical; use future_lapply for permutations (requires future.apply).

- future_seed:

  Optional seed control for future.apply (passed to future_lapply).

- progress:

  Logical; report progress via progressr (works with parallel futures).

## Value

A list with a component-level results table and summary stats.

## Details

For two-group inference, the function fits a linear model for each
component (group plus optional covariates) and uses label permutations
to build a null distribution of the component-level t-statistics. For
one-group inference, the user supplies a null distribution of \\W\\
(e.g., from permuted fits) and the observed component means are compared
to that null.

Use this when you want to identify which components show reliable group
differences in expression, rather than testing individual voxels.

- `stat` is the component-level test statistic (t-statistic for
  two-group tests; mean loading for one-group tests).

- `p_unc` is the uncorrected permutation p-value.

- `p_fwer` applies maxT family-wise error correction across components.

- `p_global` tests whether any component is significant.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires completed spatial NMF fit
  W <- matrix(rnorm(20*5), 20, 5)
  groups <- factor(rep(c("A","B"), each=10))
  result <- spatial_nmf_component_test(W=W, groups=groups, nperm=100)
} # }
```
