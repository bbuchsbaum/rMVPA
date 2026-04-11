# Get the Number of Folds

An S3 generic method to retrieve the number of folds from a
cross-validation specification object.

## Usage

``` r
get_nfolds(obj, ...)
```

## Arguments

- obj:

  A cross-validation specification object (e.g., inheriting from
  \`cross_validation\`).

- ...:

  Additional arguments passed to methods.

## Value

An integer representing the number of folds.

## Examples

``` r
cval <- kfold_cross_validation(len = 20, nfolds = 4)
get_nfolds(cval)
#> [1] 4
```
