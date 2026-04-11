# Create a Regression Result Object

This function creates a regression result object containing the observed
and predicted values, test design, test indices, and predictor.

## Usage

``` r
regression_result(
  observed,
  predicted,
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

- testind:

  A vector of indices for the test data (optional).

- test_design:

  The test design (optional).

- predictor:

  The predictor used in the regression model (optional).

## Value

A list with class attributes "regression_result",
"classification_result", and "list" containing the observed and
predicted values, test design, test indices, and predictor.

## See also

Other classification_result:
[`binary_classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/binary_classification_result.md),
[`classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/classification_result.md),
[`multiway_classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/multiway_classification_result.md)
