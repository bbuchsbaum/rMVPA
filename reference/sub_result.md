# Extract Row-wise Subset of a Result

Extract a subset of rows from a classification/regression result object.

## Usage

``` r
sub_result(x, indices)
```

## Arguments

- x:

  The input result object.

- indices:

  Row indices to extract.

## Value

A new result object with the specified rows.

## Examples

``` r
cres <- binary_classification_result(
  observed = factor(c("a","b","a")),
  predicted = factor(c("a","b","b")),
  probs = matrix(c(.8,.2,.4,.3,.7,.6), ncol=2, dimnames=list(NULL,c("a","b")))
)
sub_result(cres, 1:2)
#> 
#>  Binary Classification Result 
#> 
#> - Data Summary 
#>   - Observations:  2 
#>   - Classes:  a, b 
#>   - Test Indices:  None 
#> - Performance Metrics 
#>   - Accuracy:  1.0000 
#>   - Sensitivity:  1.0000 
#>   - Specificity:  1.0000 
#> 
```
