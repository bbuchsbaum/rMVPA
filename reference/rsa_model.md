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
  pattern_center = c("none", "stimulus_mean")
)
```

## Arguments

- dataset:

  An instance of an `mvpa_dataset`.

- design:

  An instance of an `rsa_design` created by
  [`rsa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md).

- distmethod:

  A character string specifying the method used to compute distances
  between observations. One of: `"pearson"` or `"spearman"` (defaults to
  "spearman").

- regtype:

  A character string specifying the analysis method. One of:
  `"pearson"`, `"spearman"`, `"lm"`, or `"rfit"` (defaults to
  "pearson").

- check_collinearity:

  Logical indicating whether to check for collinearity in the design
  matrix. Only applies when `regtype="lm"`. Default is TRUE.

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

## Value

An object of class `"rsa_model"` (and `"list"`), containing:

- `dataset` : the input dataset

- `design` : the RSA design

- `distmethod` : the distance method used

- `regtype` : the regression type

- `nneg` : a named list of constrained variables, if any

- `semipartial`: whether to compute semi-partial correlations

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

# Create an RSA model enforcing non-negativity for dismat2 only:
# Requires the 'glmnet' package to be installed
# rsa_mod_nneg <- rsa_model(mvpa_data, rdes, regtype="lm",
#                          nneg = list(dismat2 = TRUE))

# Create an RSA model using 'lm' but returning semi-partial correlations:
rsa_mod_sp <- rsa_model(mvpa_data, rdes, regtype = "lm",
                        semipartial = TRUE)
#> Checking design matrix for collinearity...
#> Collinearity check passed.

# Train the model using a trial-by-feature matrix
fit_params <- train_model(rsa_mod_sp, data_mat, y = NULL, indices = NULL)
# 'fit_params' = named vector of semi-partial correlations for each predictor
```
