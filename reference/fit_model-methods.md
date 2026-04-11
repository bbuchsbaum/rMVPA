# Fit Model

Fit a classification or regression model.

This method fits a multivariate pattern analysis (MVPA) model to the
provided training data.

## Usage

``` r
fit_model(
  obj,
  roi_x,
  y,
  wts,
  param,
  lev = NULL,
  last = FALSE,
  classProbs = FALSE,
  ...
)

# S3 method for class 'mvpa_model'
fit_model(
  obj,
  roi_x,
  y,
  wts,
  param,
  lev = NULL,
  last = FALSE,
  classProbs = FALSE,
  ...
)
```

## Arguments

- obj:

  An object of class `mvpa_model`.

- roi_x:

  An ROI containing the training data.

- y:

  The response vector.

- wts:

  A set of case weights.

- param:

  Tuning parameters.

- lev:

  Factor levels (for classification).

- last:

  Logical indicating if this is the last iteration.

- classProbs:

  Logical indicating if class probabilities should be returned.

- ...:

  Additional arguments to be passed to the method-specific function.

## Value

A fitted model object with additional attributes `"obsLevels"` and
`"problemType"`.

## Examples

``` r
# \donttest{
  if (requireNamespace("sda", quietly = TRUE)) {
    ds <- gen_sample_dataset(
      D = c(6, 6, 6), nobs = 20,
      response_type = "categorical",
      data_mode = "image", nlevels = 2
    )
    mdl <- load_model("sda_notune")
    mspec <- mvpa_model(
      model = mdl,
      dataset = ds$dataset,
      design  = ds$design,
      model_type = "classification"
    )
    grid <- tune_grid(mspec, ds$dataset$train_data, ds$design$y_train, len = 1)
    fit  <- fit_model(mspec, ds$dataset$train_data,
                     ds$design$y_train, wts = NULL, param = grid)
  }
# }
```
