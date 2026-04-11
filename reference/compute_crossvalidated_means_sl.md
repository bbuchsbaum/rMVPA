# Compute Cross-Validated Condition Means within a Searchlight

This helper function calculates the mean activation pattern for each
condition using data from other cross-validation folds.

## Usage

``` r
compute_crossvalidated_means_sl(
  sl_data,
  mvpa_design,
  cv_spec,
  estimation_method = "average",
  whitening_matrix_W = NULL,
  return_folds = FALSE
)
```

## Arguments

- sl_data:

  A numeric matrix (samples x voxels/vertices) containing the data for
  the current searchlight.

- mvpa_design:

  The `mvpa_design` object associated with the dataset, containing
  condition labels and block information.

- cv_spec:

  An object describing the cross-validation scheme, typically created by
  functions like `\link{blocked_cross_validation}`,
  `\link{twofold_blocked_cross_validation}`,
  `\link{kfold_cross_validation}`, etc. (inheriting from
  `cross_validation`). This object determines how training/test folds
  are defined.

- estimation_method:

  Character string specifying the method to estimate means. Currently
  supported: `"average"` (simple mean of training samples per
  condition).

  - `"average"`: Simple mean of training samples per condition.

  - `"L2_norm"`: Identical to `"average"` but each condition pattern
    (row) is finally scaled to unit L2 norm. Useful when you need to
    equalise overall pattern energy across conditions before RSA.

  - `"crossnobis"`: Applies a pre-computed whitening matrix (see
    \`whitening_matrix_W\`) to the average pattern of each condition
    within each cross-validation training fold, before averaging these
    whitened patterns across folds. This aims to produce
    noise-normalized condition representations.

  Default is `"average"`.

- whitening_matrix_W:

  Optional V x V numeric matrix, where V is the number of
  voxels/features in \`sl_data\`. This matrix should be the whitening
  transformation (e.g., Sigma_noise^(-1/2)) derived from GLM residuals.
  Required and used only if \`estimation_method = "crossnobis"\`.

- return_folds:

  Logical, if TRUE, the function returns a list containing both the
  overall mean estimate (\`mean_estimate\`) and an array of per-fold
  estimates (\`fold_estimates\`). If FALSE (default), only the overall
  mean estimate is returned.

## Value

If \`return_folds = FALSE\` (default): A numeric matrix (K x V_sl) where
K is the number of conditions and V_sl is the number of voxels/vertices
in the searchlight. Each row represents the cross-validated mean pattern
for condition k. If \`return_folds = TRUE\`: A list with two elements:

- \`mean_estimate\`:

  The K x V_sl matrix described above.

- \`fold_estimates\`:

  A K x V_sl x M array, where M is the number of folds, containing the
  mean estimate for each condition from each fold.

## Examples

``` r
if (FALSE) { # \dontrun{
  # See vignette for cross-validated mean computation
} # }
```
