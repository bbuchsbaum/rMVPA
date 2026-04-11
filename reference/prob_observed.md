# Probability of Observed Class

Extract the predicted probability for the observed class.

## Usage

``` r
prob_observed(x)
```

## Arguments

- x:

  The object from which to extract the probability.

## Value

A vector of predicted probabilities.

## Examples

``` r
cres <- binary_classification_result(
  observed = factor(c("a","b")),
  predicted = factor(c("a","b")),
  probs = matrix(c(.8,.2,.3,.7), ncol=2, dimnames=list(NULL,c("a","b")))
)
prob_observed(cres)
#> [1] 0.8 0.7
```
