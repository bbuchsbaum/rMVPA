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

[`regression_result`](http://bbuchsbaum.github.io/rMVPA/reference/regression_result.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # See performance() generic for examples
} # }
```
