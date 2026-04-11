# Classification results for binary outcome

Constructs a binary classification result object based on the observed
and predicted values, as well as other optional parameters.

## Usage

``` r
binary_classification_result(
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

A binary classification result object, with the class attribute set to
"binary_classification_result".

## See also

Other classification_result:
[`classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/classification_result.md),
[`multiway_classification_result()`](http://bbuchsbaum.github.io/rMVPA/reference/multiway_classification_result.md),
[`regression_result()`](http://bbuchsbaum.github.io/rMVPA/reference/regression_result.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  obs <- factor(c("a", "b", "a", "b"))
  pred <- factor(c("a", "b", "b", "a"))
  probs <- matrix(c(0.8, 0.2, 0.7, 0.3, 0.4, 0.6, 0.3, 0.7),
                  ncol=2, dimnames=list(NULL, c("a","b")))
  res <- binary_classification_result(obs, pred, probs)
} # }
```
