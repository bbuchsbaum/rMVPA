# Format Result Object

Format Result Object

## Usage

``` r
format_result(obj, result, error_message, ...)

# S3 method for class 'mvpa_model'
format_result(obj, result, error_message = NULL, context, ...)

# S3 method for class 'feature_rsa_model'
format_result(obj, result, error_message = NULL, context, ...)
```

## Arguments

- obj:

  The result object to be formatted.

- result:

  The result object to be formatted.

- error_message:

  An optional error message.

- ...:

  Additional arguments to be passed to the method-specific function.

- context:

  Optional contextual metadata (e.g., ROI details) supplied to method
  implementations.

## Value

A formatted result object, typically a tibble row.

## Examples

``` r
# \donttest{
dataset <- gen_sample_dataset(D = c(6, 6, 6), nobs = 20,
                              response_type = "categorical",
                              data_mode = "image")
cval <- blocked_cross_validation(dataset$design$block_var)
model <- load_model("sda_notune")
mspec <- mvpa_model(model, dataset$dataset, dataset$design,
                    "classification", crossval = cval)
# Typically called internally during processing
format_result(mspec, result = NULL, error_message = "example",
              context = list(test = 1, ytest = factor("a")))
#> # A tibble: 1 × 6
#>   class  probs  y_true    fit    error error_message
#>   <list> <list> <list>    <list> <lgl> <chr>        
#> 1 <NULL> <NULL> <fct [1]> <NULL> TRUE  example      
# }
```
