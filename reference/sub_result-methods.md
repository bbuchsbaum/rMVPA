# Subset Multiway Classification Result

This function subsets a multiway classification result based on the
provided indices.

This function subsets a binary classification result based on the
provided indices.

## Usage

``` r
# S3 method for class 'multiway_classification_result'
sub_result(x, indices)

# S3 method for class 'binary_classification_result'
sub_result(x, indices)
```

## Arguments

- x:

  An object of class `binary_classification_result` containing the
  binary classification results.

- indices:

  The set of indices used to subset the results.

## Value

A `multiway_classification_result` object containing the subset of
results specified by the indices.

A `binary_classification_result` object containing the subset of results
specified by the indices.

## Examples

``` r
if (FALSE) { # \dontrun{
  # S3 method - called via sub_result generic
  cres <- multiway_classification_result(
    factor(c("a","b","c")), factor(c("a","b","c")),
    matrix(runif(9), 3, 3, dimnames=list(NULL, c("a","b","c")))
  )
  subset <- sub_result(cres, 1:2)
} # }
```
