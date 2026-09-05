# Construct an RSA (Representational Similarity Analysis) model

This function creates an RSA model object by taking an MVPA
(Multi-Variate Pattern Analysis) dataset and an RSA design.

## Usage

``` r
rsa_model(
  dataset,
  design,
  distmethod = "spearman",
  regtype = "pearson",
  check_collinearity = TRUE,
  nneg = NULL,
  semipartial = FALSE,
  pattern_center = c("none", "stimulus_mean"),
  return_fingerprint = FALSE,
  fingerprint_method = c("pearson", "spearman"),
  fingerprint_basis = c("pca", "qr")
)
```

## Arguments

- dataset:

  An instance of an `mvpa_dataset`.

- design:

  An instance of an `rsa_design` created by
  [`rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md).

- distmethod:

  A character string specifying the method used to compute distances
  between observations. One of: `"pearson"` or `"spearman"` (defaults to
  "spearman").

- regtype:

  A character string specifying the analysis method. One of:
  `"pearson"`, `"spearman"`, `"lm"`, or `"rfit"` (defaults to
  "pearson").

- check_collinearity:

  Logical. When `regtype = "lm"`, stop if two predictor RDMs correlate
  above 0.99 or the design matrix is rank deficient. When `regtype` is
  `"lm"` or `"rfit"`, also warn if a model predictor has no usable
  variation or fewer than 10 effective items (see the section on
  effective support). Default is TRUE.

- nneg:

  A named list of variables (predictors) for which non-negative
  regression coefficients should be enforced (only if `regtype="lm"`).
  Defaults to `NULL` (no constraints).

- semipartial:

  Logical indicating whether to compute semi-partial correlations in the
  `"lm"` case (only if `nneg` is not used). Defaults to `FALSE`.

- pattern_center:

  Optional pattern-centering method applied to the stimulus-by-voxel
  matrix before distances are computed. Use `"stimulus_mean"` to
  subtract the across-stimulus mean pattern (Hanson-style). Default is
  `"none"`.

- return_fingerprint:

  Logical; if `TRUE`, project the per-unit neural pair-dissimilarity
  vector onto an orthonormal basis of the signal model RDM subspace and
  return the standardized model-space score vector alongside the
  standard outputs. For
  [`pair_rsa_design`](https://bbuchsbaum.github.io/rMVPA/reference/pair_rsa_design.md)
  objects, predictors supplied through `nuisance` are excluded from this
  fingerprint basis. The basis is cached on the model spec and reused
  across ROIs/searchlight units. Default `FALSE`.

- fingerprint_method:

  Standardization method used when `return_fingerprint = TRUE`. One of
  `"pearson"` (centered/scaled) or `"spearman"` (rank-then-standardize).
  Default `"pearson"`.

- fingerprint_basis:

  Basis used to span the model RDM subspace when
  `return_fingerprint = TRUE`: `"pca"` (default) or `"qr"`.

## Value

An object of class `"rsa_model"` (and `"list"`), containing:

- `dataset` : the input dataset

- `design` : the RSA design

- `distmethod` : the distance method used

- `regtype` : the regression type

- `nneg` : a named list of constrained variables, if any

- `semipartial`: whether to compute semi-partial correlations

- `design_diagnostics`: an `rsa_design_diagnostics` object (see
  [`rsa_design_diagnostics`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_design_diagnostics.md))

## Effective support of the design

RDM entries are not independent observations. Every item enters
`n_items - 1` pairs, so the `n_pairs` entries of a vectorised RDM carry
on the order of `n_items` independent pieces of information, and a
regression coefficient's variance is inflated further by the predictor's
collinearity with the others (its VIF). The constructor computes
`n_items / VIF` for every predictor and stores the result as
`design_diagnostics` (see
[`rsa_design_diagnostics`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_design_diagnostics.md)).
When `regtype` is `"lm"` or `"rfit"` and a model predictor falls below
10 effective items, it warns: that coefficient will vary across ROIs
largely through noise. This is a heuristic screen, not a test.

For inference on RSA maps use
[`run_permutation_searchlight`](https://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md),
which permutes item labels and so carries the same dependence into the
null. For regression with multiple design predictors (including nuisance
predictors), this supports only an explicit joint no-association null
via `permutation_control(rsa_null = "joint")`; it does not provide
conditional inference on individual coefficients. The per-ROI t-values
returned by `regtype = "lm"` use pair-based degrees of freedom and are
anti-conservative.

## Examples

``` r
# Create a random MVPA dataset (image data)
arr  <- array(rnorm(100 * 5), c(5, 5, 4, 5))   # 5 voxels x 5 voxels x 4 slices x 5 observations
sp   <- neuroim2::NeuroSpace(c(5, 5, 4, 5))
vec  <- neuroim2::NeuroVec(arr, sp)
mask <- neuroim2::LogicalNeuroVol(array(1, c(5, 5, 4)), neuroim2::NeuroSpace(c(5, 5, 4)))
mvpa_data <- mvpa_dataset(train_data = vec, mask = mask)

# Create two random RDMs (distance matrices) over the 5 observations
data_mat  <- matrix(rnorm(5 * 10), 5, 10)
dismat1   <- dist(data_mat)
dismat2   <- dist(matrix(rnorm(5 * 10), 5, 10))
rdes <- rsa_design(~ dismat1 + dismat2,
                   list(dismat1 = dismat1, dismat2 = dismat2))

# Create an RSA model with standard 'lm' (returns t-values):
rsa_mod <- rsa_model(mvpa_data, rdes, regtype = "lm")
#> Checking design matrix for collinearity...
#> Collinearity check passed.
#> Warning: RSA design has only 5 items, which is ~4.1 effective items for predictor 'dismat1' (VIF = 1.2). Its coefficient will be unstable across ROIs. Add items or pool runs. RDM entries share items, so the 10 pairs carry roughly 5 independent observations. Per-ROI t-values overstate the evidence. See run_permutation_searchlight() for the supported RSA null hypotheses; conditional coefficient inference with other predictors is not implemented.

# Create an RSA model enforcing non-negativity for dismat2 only:
# Requires the 'glmnet' package to be installed
# rsa_mod_nneg <- rsa_model(mvpa_data, rdes, regtype="lm",
#                          nneg = list(dismat2 = TRUE))

# Create an RSA model using 'lm' but returning semi-partial correlations:
rsa_mod_sp <- rsa_model(mvpa_data, rdes, regtype = "lm",
                        semipartial = TRUE)
#> Checking design matrix for collinearity...
#> Collinearity check passed.
#> Warning: RSA design has only 5 items, which is ~4.1 effective items for predictor 'dismat1' (VIF = 1.2). Its coefficient will be unstable across ROIs. Add items or pool runs. RDM entries share items, so the 10 pairs carry roughly 5 independent observations. Per-ROI t-values overstate the evidence. See run_permutation_searchlight() for the supported RSA null hypotheses; conditional coefficient inference with other predictors is not implemented.

# Train the model using a trial-by-feature matrix
fit_params <- train_model(rsa_mod_sp, data_mat, y = NULL, indices = NULL)
# 'fit_params' = named vector of semi-partial correlations for each predictor
```
