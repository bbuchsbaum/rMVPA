# Calculate the Predicted Class from Probability Matrix

This function calculates the predicted class from a matrix of predicted
probabilities. The class with the highest probability is selected as the
predicted class.

## Usage

``` r
predicted_class(prob)
```

## Arguments

- prob:

  A matrix of predicted probabilities with column names indicating the
  classes.

## Value

A vector of predicted classes corresponding to the highest probability
for each row in the input matrix.

## Examples

``` r
prob <- matrix(c(0.2, 0.8,
                 0.6, 0.4),
               nrow = 2, byrow = TRUE,
               dimnames = list(NULL, c("A", "B")))
predicted_class(prob)
```
