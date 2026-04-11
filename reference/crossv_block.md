# Block Cross-Validation Data Preparation

This function prepares the data for block cross-validation by dividing
the dataset based on the provided block variable. It creates subsets of
training and testing data for each block without performing any analysis
or fitting models.

## Usage

``` r
crossv_block(data, y, block_var, id = ".id")
```

## Arguments

- data:

  A data frame containing the training data.

- y:

  A response vector.

- block_var:

  An integer vector defining the cross-validation blocks.

- id:

  A character string specifying the identifier for the output data
  frame.

## Value

A tibble containing the training and testing data, response vectors, and
block IDs for each fold.

## Examples

``` r
X <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
y <- rep(letters[1:4], 25)
block_var <- rep(1:4, each = 25)
cv <- crossv_block(X, y, block_var)
```
