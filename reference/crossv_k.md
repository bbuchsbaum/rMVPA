# K-fold Cross-Validation Data Preparation

This function prepares the data for k-fold cross-validation by dividing
the dataset into k folds. It creates subsets of training and testing
data for each fold without performing any analysis or fitting models.

## Usage

``` r
crossv_k(data, y, k = 5, id = ".id")
```

## Arguments

- data:

  A data frame containing the training data.

- y:

  A response vector.

- k:

  An integer specifying the number of folds for cross-validation.

- id:

  A character string specifying the identifier for the output data
  frame.

## Value

A tibble containing the training and testing data, response vectors, and
fold IDs for each fold.

## Examples

``` r
data <- iris[,-5]
y <- iris$Species
result <- crossv_k(data, y, k = 5)
```
