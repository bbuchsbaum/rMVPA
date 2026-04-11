# Compute Performance for an Object

Delegates calculation of performance metrics to the appropriate method.

## Usage

``` r
compute_performance(obj, result)

# S3 method for class 'mvpa_model'
compute_performance(obj, result)
```

## Arguments

- obj:

  Model specification or object capable of computing performance.

- result:

  The classification/regression result to evaluate.

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
dummy <- list(performance = performance)
class(dummy) <- "mvpa_model"
compute_performance(dummy, cres)
#> Accuracy      AUC 
#>        1        1 
```
