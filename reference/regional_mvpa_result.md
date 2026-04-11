# Create a `regional_mvpa_result` instance

Constructs a regional MVPA result object that stores the results of MVPA
analysis in a specific region.

## Usage

``` r
regional_mvpa_result(
  model_spec,
  performance_table,
  prediction_table,
  vol_results,
  fits = fits,
  pooled_prediction_table = NULL,
  pooled_performance = NULL
)
```

## Arguments

- model_spec:

  A model specification object.

- performance_table:

  A data frame with performance measures.

- prediction_table:

  A data frame with prediction results.

- vol_results:

  A list of voxel-level results.

- fits:

  Optional model fits.

- pooled_prediction_table:

  Optional pooled prediction table (e.g., ROI-averaged or stacked).

- pooled_performance:

  Optional named numeric vector of pooled performance metrics.

## Value

A `regional_mvpa_result` object.

## Examples

``` r
# Create example inputs
model_spec <- list(dataset = "Example dataset")
performance_table <- data.frame(accuracy = c(0.8, 0.85))
prediction_table <- data.frame(observed = factor(rep(letters[1:2], 5)),
                                predicted = factor(rep(letters[1:2], 5)))
vol_results <- list(vol1 = "Example vol_result 1", vol2 = "Example vol_result 2")
fits <- list(fit1 = "Example fit 1", fit2 = "Example fit 2")

# Construct a regional_mvpa_result
regional_result <- regional_mvpa_result(model_spec, performance_table,
                                        prediction_table, vol_results, fits = fits)
```
