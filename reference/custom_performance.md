# Apply Custom Performance Metric to Prediction Result

This function applies a user-supplied performance metric to a prediction
result object. The custom function should take a single result object
(e.g., a `classification_result`) and return a named numeric vector of
scalar metrics.

## Usage

``` r
custom_performance(x, custom_fun, split_list = NULL)
```

## Arguments

- x:

  The prediction result object.

- custom_fun:

  The function used to compute performance metrics, i.e.,
  `custom_fun(x)`; must return a numeric vector (named if multiple
  metrics are returned).

- split_list:

  An optional named list of splitting groups. If provided, metrics are
  computed on the full result and separately within each group, with
  group-specific metrics suffixed by `"_<group>"`.

## Value

A named numeric vector with the calculated custom performance metric(s).

## Details

When a design includes `split_by`, the corresponding split groups
(typically `mvpa_design$split_groups`) can be passed as `split_list`. If
the result object contains `testind`, indices in `split_list` are
interpreted in the design space and mapped to rows of the result.

## Examples

``` r
cres <- binary_classification_result(
  observed  = factor(c("A", "B")),
  predicted = factor(c("A", "A")),
  probs = matrix(c(0.9, 0.1,
                   0.6, 0.4),
                 ncol = 2, byrow = TRUE,
                 dimnames = list(NULL, c("A", "B")))
)
acc_fun <- function(x) c(acc = mean(x$observed == x$predicted))
custom_performance(cres, acc_fun)
#> acc 
#> 0.5 
```
