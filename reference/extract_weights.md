# Extract Raw Model Weights

Extract the raw weight matrix from a fitted model object. Returns a P x
D numeric matrix (features x discriminant directions).

## Usage

``` r
extract_weights(object, ...)

# S3 method for class 'sda'
extract_weights(object, ...)

# S3 method for class 'glmnet'
extract_weights(object, ...)

# S3 method for class 'spacenet_fit'
extract_weights(object, ...)

# S3 method for class 'model_fit'
extract_weights(object, ...)

# Default S3 method
extract_weights(object, ...)
```

## Arguments

- object:

  A fitted model object (e.g., from `sda`, `glmnet`).

- ...:

  Additional arguments passed to methods.

## Value

A numeric matrix of dimension P x D.

## Examples

``` r
# \donttest{
  if (requireNamespace("sda", quietly = TRUE)) {
    ds <- gen_sample_dataset(c(5,5,5), 20, nlevels=2)
    mdl <- load_model("sda_notune")
    mspec <- mvpa_model(mdl, ds$dataset, ds$design, "classification")
    vox <- which(ds$dataset$mask > 0)
    X <- neuroim2::series(ds$dataset$train_data, vox)
    fit <- train_model(mspec, X, ds$design$y_train, indices=vox)
    w <- extract_weights(fit)
  }
# }
```
