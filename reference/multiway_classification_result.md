# Create a Multiway Classification Result Object

This function creates a multiway classification result object containing
the observed and predicted values, class probabilities, test design,
test indices, and predictor.

## Usage

``` r
multiway_classification_result(
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

  A vector of observed values.

- predicted:

  A vector of predicted values.

- probs:

  A matrix of class probabilities.

- testind:

  A vector of indices for the test data (optional).

- test_design:

  The test design (optional).

- predictor:

  The predictor used in the multiway classification model (optional).

## Value

A list with class attributes "multiway_classification_result",
"classification_result", and "list" containing the observed and
predicted values, class probabilities, test design, test indices, and
predictor.

## See also

Other classification_result:
[`binary_classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/binary_classification_result.md),
[`classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/classification_result.md),
[`regression_result()`](http://bbuchsbaum.github.io/rMVPA/reference/regression_result.md)
