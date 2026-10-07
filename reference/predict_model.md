# Predict Model Output

Generic function to predict outcomes from a fitted model object using
new data.

## Usage

``` r
predict_model(object, fit, newdata, ...)
```

## Arguments

- object:

  A fitted model object for which a prediction method is defined.

- fit:

  The fitted model object, often returned by \`train_model\`. (Note: For
  some models, \`object\` itself might be the fit).

- newdata:

  New data (e.g., a matrix or data frame) for which to make predictions.
  The structure should be compatible with what the model was trained on.

- ...:

  Additional arguments passed to specific prediction methods.

## Value

Predictions whose structure depends on the specific method (e.g., a
vector, matrix, or data frame).

## Examples

``` r
set.seed(1)
n <- 24
feats <- matrix(rnorm(n * 4), n, 4)    # stimulus features
brain <- matrix(rnorm(n * 10), n, 10)  # ROI patterns (trials x voxels)
dset <- gen_sample_dataset(c(3, 3, 3), n, blocks = 3)
fdes <- feature_rsa_design(F = feats, labels = seq_len(n), max_comps = 3)
mspec <- feature_rsa_model(dset$dataset, fdes, method = "pls",
  crossval = blocked_cross_validation(dset$design$block_var))
fit <- train_model(mspec, brain, feats, indices = seq_len(ncol(brain)))
preds <- predict_model(mspec, fit, feats[1:5, ])
dim(preds)
#> [1]  5 10
```
