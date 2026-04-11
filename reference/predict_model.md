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
if (FALSE) { # \dontrun{
  preds <- predict_model(model_spec, fitted_model, new_data)
} # }
```
