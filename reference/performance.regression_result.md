# Calculate Performance Metrics for Regression Result

This function calculates performance metrics for a regression result
object, including R-squared, Root Mean Squared Error (RMSE), and
Spearman correlation.

## Usage

``` r
# S3 method for class 'regression_result'
performance(x, split_list = NULL, ...)
```

## Arguments

- x:

  A `regression_result` object.

- split_list:

  Optional named list of split index groups for computing metrics on
  sub-groups in addition to the full result.

- ...:

  extra args (not used).

## Value

A named vector with the calculated performance metrics: R-squared, RMSE,
and Spearman correlation.

## Details

The function calculates the following performance metrics for the given
regression result object: - R-squared: proportion of variance in the
observed data that is predictable from the fitted model. - RMSE: root
mean squared error, a measure of the differences between predicted and
observed values. - Spearman correlation: a measure of the monotonic
relationship between predicted and observed values.

## See also

[`regression_result`](https://bbuchsbaum.github.io/rMVPA/reference/regression_result.md)

## Examples

``` r
res <- regression_result(observed = c(1.0, 2.1, 2.9, 4.2, 5.1),
                         predicted = c(1.2, 1.9, 3.1, 3.8, 5.0),
                         testind = 1:5)
performance(res)
#>        R2      RMSE  spearcor 
#> 0.9786905 0.2408319 1.0000000 
```
