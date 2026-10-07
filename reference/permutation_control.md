# Create a Permutation Control Object

Specifies all tuning parameters for
[`run_permutation_searchlight`](https://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md).

## Usage

``` r
permutation_control(
  n_perm = 5L,
  shuffle = c("within_block", "circular_shift", "global"),
  null_method = c("adjusted", "global"),
  adjust_by = c("nfeatures", "redundancy", "both"),
  n_bins = 5L,
  subsample = 0.1,
  stratify_subsample = TRUE,
  correction = c("fdr", "none"),
  diagnose = TRUE,
  seed = NULL,
  perm_strategy = c("iterate", "searchlight"),
  rsa_null = c("individual", "joint")
)
```

## Arguments

- n_perm:

  Integer \>= 1. Number of permutations to run.

- shuffle:

  Character. How to permute labels: `"within_block"` (default) shuffles
  within each block, `"circular_shift"` shifts the label sequence within
  each block, `"global"` shuffles all labels ignoring block structure.
  Circular shifts include zero in each block; leaving some or all blocks
  unchanged is a valid null draw.

- null_method:

  Character. How to build the null distribution: `"adjusted"` conditions
  on covariate bins (default), `"global"` uses one global null.

- adjust_by:

  Character. Which covariate(s) to condition on: `"nfeatures"`
  (default), `"redundancy"`, or `"both"`.

- n_bins:

  Integer \>= 2. Number of quantile bins for covariate stratification.

- subsample:

  Fraction (0, 1\] or integer count of searchlight centers to use for
  permutation runs. Only used when `perm_strategy = "iterate"` (see
  below).

- stratify_subsample:

  Logical. If `TRUE` (default), subsample centers proportionally from
  nfeatures quantile bins.

- correction:

  Character. Multiple-comparison correction: `"fdr"`
  (Benjamini-Hochberg, default) or `"none"`.

- diagnose:

  Logical. If `TRUE` (default), run null diagnostics.

- seed:

  Optional integer random seed.

- perm_strategy:

  Character. Controls how each permutation pass is executed. Eligible
  exact engines reuse preparation across permutations for built-in
  blocked CV or explicit custom splits. Other CV specifications retain
  their per-permutation fold generation.

  `"iterate"` (default)

  :   Each permutation evaluates a **subsampled** set of centers using a
      prepared engine when eligible, or
      [`mvpa_iterate`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
      otherwise.

      **When to use**: slow classifiers, large brains, limited compute.
      The `subsample` parameter controls how many centers are evaluated
      per permutation.

      Null pool size: `n_perm * n_subsampled_centers`.

  `"searchlight"`

  :   Each permutation evaluates the **full brain** using a prepared
      engine when eligible, or
      [`run_searchlight`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
      otherwise. The latter retains standard engine dispatch and
      user-defined `run_searchlight.<class>` methods.

      **When to use**: models with a fast searchlight engine, or when
      you want the richest possible null distribution. Since the full
      brain is computed anyway, *all* centers contribute to the null
      (the `subsample` parameter is ignored and a note is logged).

      Null pool size: `n_perm * all_centers`.

- rsa_null:

  Null hypothesis for `rsa_model`: `"individual"` (default) tests
  marginal associations for correlation models, or the coefficient of a
  regression with only one design predictor. Regression models with
  multiple predictors (including nuisance predictors) require `"joint"`:
  no association with *any* design predictor. Raw item permutations
  cannot test an individual regression coefficient conditional on the
  others. Ignored for other model classes.

## Value

An S3 object of class `"permutation_control"`.

## Examples

``` r
# Default: subsampled iterator (safe for any model)
pc <- permutation_control(n_perm = 100, shuffle = "within_block",
                          subsample = 0.1, seed = 42L)

# Full-brain strategy (benefits from fast engines, richer null)
pc2 <- permutation_control(n_perm = 20, perm_strategy = "searchlight",
                           seed = 42L)
print(pc2)
#> Permutation Control Settings
#>   perm_strategy     : searchlight 
#>   n_perm            : 20 
#>   shuffle           : within_block 
#>   rsa_null          : individual 
#>   null_method       : adjusted 
#>   adjust_by         : nfeatures 
#>   n_bins            : 5 
#>   subsample         : (ignored - full brain per perm)
#>   correction        : fdr 
#>   diagnose          : TRUE 
#>   seed              : 42 
```
