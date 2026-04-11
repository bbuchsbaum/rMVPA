# Create a `classification_result` instance

Constructs a classification result object based on the observed and
predicted values, as well as other optional parameters.

## Usage

``` r
classification_result(
  observed,
  predicted,
  probs,
  testind = NULL,
  test_design = NULL,
  predictor = NULL
)
```

## Arguments

- observed:

  A vector of observed or true values.

- predicted:

  A vector of predicted values.

- probs:

  A `matrix` of predicted probabilities, with one column per level.

- testind:

  The row indices of the test observations (optional).

- test_design:

  An optional design for the test data.

- predictor:

  An optional predictor object.

## Value

A classification result object, which can be one of:
`regression_result`, `binary_classification_result`, or
`multiway_classification_result`.

## See also

Other classification_result:
[`binary_classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/binary_classification_result.md),
[`multiway_classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/multiway_classification_result.md),
[`regression_result()`](http://bbuchsbaum.github.io/rMVPA/reference/regression_result.md)

## Examples

``` r
# A vector of observed values
yobs <- factor(rep(letters[1:4], 5))

# Predicted probabilities
probs <- data.frame(a = runif(1:20), b = runif(1:20), c = runif(1:20), d = runif(1:20))
probs <- sweep(probs, 1, rowSums(probs), "/")

# Get the max probability per row and use this to determine the predicted class
maxcol <- max.col(probs)
predicted <- levels(yobs)[maxcol]

# Construct a classification result
cres <- classification_result(yobs, predicted, probs)

# Compute default performance measures (Accuracy, AUC)
performance(cres)
#>    Accuracy         AUC 
#>  0.25000000 -0.04666667 
```
