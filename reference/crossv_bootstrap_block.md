# Block Bootstrap Cross-Validation Data Preparation

This function prepares the data for block bootstrap cross-validation by
dividing the dataset based on the provided block variable. It creates
subsets of training and testing data for each block using bootstrap
sampling within the training blocks, without performing any analysis or
fitting models.

## Usage

``` r
crossv_bootstrap_block(
  data,
  y,
  block_var,
  nreps = 5,
  id = ".id",
  weights = NULL
)
```

## Arguments

- data:

  A data frame containing the training data.

- y:

  A response vector.

- block_var:

  An integer vector defining the cross-validation blocks.

- nreps:

  An integer specifying the number of bootstrap repetitions.

- id:

  A character string specifying the identifier for the output data
  frame.

- weights:

  An optional numeric vector of weights to be used for bootstrap
  sampling.

## Value

A tibble containing the training and testing data, response vectors, and
block IDs for each fold.

## Details

The function first checks if the length of the \`block_var\` vector
matches the length of the response vector \`y\`. It then creates a list
of block indices and ensures there is more than one block to bootstrap.
If weights are provided, the function splits the weights according to
the block variable.

The function performs bootstrap sampling within the training blocks but
keeps the test set fixed. For each block, it generates a list of
training indices using bootstrap sampling and creates the corresponding
training and testing data sets.
