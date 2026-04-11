# Compute Performance Metrics

Generic function to compute performance metrics from result objects.

## Usage

``` r
performance(x, ...)
```

## Arguments

- x:

  Result object from a classification or regression analysis.

- ...:

  Additional arguments passed to methods.

## Value

Named numeric vector of performance metrics.

## Examples

``` r
cres <- binary_classification_result(
  observed  = factor(c("a", "b")),
  predicted = factor(c("a", "b")),
  probs     = matrix(c(0.8, 0.2, 0.3, 0.7), ncol = 2,
                     dimnames = list(NULL, c("a", "b")))
)
performance(cres)
#> Accuracy      AUC 
#>        1        1 
```
