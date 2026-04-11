# Train a classification, regression, or representational model.

This is a generic function that trains a model based on the provided
model specification object. Different model types will have different
methods implemented with specific parameters.

This function trains a multivariate analysis of variance (MANOVA) model
using the specified design.

This function implements the core logic for the MS-ReVE analysis within
a single searchlight or region.

## Usage

``` r
train_model(obj, ...)

# S3 method for class 'manova_model'
train_model(obj, train_dat, y, indices, ...)

# S3 method for class 'mvpa_model'
train_model(obj, train_dat, y, indices, wts = NULL, ...)

# S3 method for class 'rsa_model'
train_model(obj, train_dat, y, indices, ...)

# S3 method for class 'vector_rsa_model'
train_model(obj, train_dat, y, indices, ...)

# S3 method for class 'feature_rsa_model'
train_model(obj, X, y, indices, ...)

# S3 method for class 'contrast_rsa_model'
train_model(obj, sl_data, sl_info, cv_spec, ...)
```

## Arguments

- obj:

  An object of class `contrast_rsa_model`.

- ...:

  Additional arguments (currently ignored).

- train_dat:

  A data frame or matrix representing the training subset (e.g., voxel
  intensities).

- y:

  Feature matrix used for RSA (samples x features).

- indices:

  Spatial indices associated with the training data.

- wts:

  Optional class weights (if the underlying model supports it).

- X:

  Brain data (samples x voxels).

- sl_data:

  The data matrix for the current searchlight (samples x voxels).

- sl_info:

  A list containing information about the current searchlight, including
  `center_local_id`.

- cv_spec:

  The cross-validation specification.

## Value

A trained model object. The exact return value depends on the specific
method implementation.

A named numeric vector of -log(p-values) for each predictor in the
MANOVA model.

A model fit object containing the trained model, its fit, the model type
(classification or regression), the best tuning parameters, the voxel
indices, and the feature mask.

Depending on `obj$regtype`:

- `"lm"` + no constraints + `obj$semipartial=TRUE`: semi-partial
  correlations

- `"lm"` + no constraints + `obj$semipartial=FALSE`: T-values of each
  predictor

- `"lm"` + `nneg` constraints: raw coefficients from constrained
  `glmnet`

- `"rfit"`: robust regression coefficients

- `"pearson"` or `"spearman"`: correlation coefficients

A structure containing "scores" or similar second-order similarity
results.

A list containing RSA metrics and, if requested, permutation results.

A named list where each element corresponds to a requested
\`output_metric\` from the \`obj\$output_metric\` vector. Each element
is:

- For metrics like "beta_delta", "beta_only", "delta_only": A Q-length
  named vector where values are indexed by contrast and names match the
  contrast_matrix column names (Q = number of contrasts)

- For metrics like "recon_score", "composite": A single numeric value

The list will have an attribute "na_reason" if any metric calculation
failed, which can be used for diagnostics.

For example, if \`obj\$output_metric = c("beta_delta", "recon_score")\`,
the returned list will have two elements: \`\$beta_delta\` (a Q-length
vector) and \`\$recon_score\` (a single value).

A named list where each element corresponds to a requested metric from
`obj$output_metric`.

## Examples

``` r
# \donttest{
  # Generate a small sample dataset for classification
  dset_info <- gen_sample_dataset(
    D = c(8, 8, 8),
    nobs = 20,
    response_type = "categorical",
    data_mode = "image",
    nlevels = 2
  )

  # Create a cross-validation specification
  cval <- blocked_cross_validation(dset_info$design$block_var)

  # Load a simple classifier
  sda_model <- load_model("sda_notune")

  # Create an MVPA model specification
  mspec <- mvpa_model(
    model = sda_model,
    dataset = dset_info$dataset,
    design = dset_info$design,
    model_type = "classification",
    crossval = cval
  )

  # Extract voxel-wise data as a numeric matrix for training
  vox <- which(dset_info$dataset$mask > 0)
  X   <- neuroim2::series(dset_info$dataset$train_data, vox)

  # Train the model on the voxel matrix
  fit <- train_model(
    mspec,
    X,
    dset_info$design$y_train,
    indices = vox
  )
# }
# This example shows the structure of the returned list but doesn't actually run the function
# For a multi-metric model: output_metric = c("beta_delta", "recon_score", "beta_only")
```
