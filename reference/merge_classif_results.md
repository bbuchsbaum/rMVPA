# Merge Multiple Classification/Regression Results

This function merges two or more classification/regression result
objects.

## Usage

``` r
merge_classif_results(x, ...)
```

## Arguments

- x:

  The first classification/regression result object.

- ...:

  Additional classification/regression result objects.

## Value

A single merged classification/regression result object.

## Examples

``` r
cres1 <- binary_classification_result(
  observed = factor(c("a","b")),
  predicted = factor(c("a","b")),
  probs = matrix(c(.8,.2,.3,.7), ncol=2, dimnames=list(NULL,c("a","b")))
)
cres2 <- binary_classification_result(
  observed = factor(c("b","a")),
  predicted = factor(c("b","a")),
  probs = matrix(c(.2,.8,.7,.3), ncol=2, dimnames=list(NULL,c("a","b")))
)
merged <- merge_classif_results(cres1, cres2)
```
