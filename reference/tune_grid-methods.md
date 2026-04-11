# Extract Tuning Grid

Returns the parameter grid used to tune a model.

## Usage

``` r
tune_grid(obj, x, y, len)
```

## Arguments

- obj:

  A model or model specification.

- x:

  Training data.

- y:

  Response variable.

- len:

  Number of parameter sets to generate.

## Value

A data frame of tuning parameter combinations.

## Examples

``` r
# \donttest{
ds  <- gen_sample_dataset(D = c(5, 5, 5), nobs = 10)
mdl <- load_model("sda_notune")
mspec <- mvpa_model(mdl, ds$dataset, ds$design, model_type = "classification")
tune_grid(mspec, ds$dataset$train_data, ds$design$y_train, len = 1)
#>   parameter
#> 1      none
# }
```
