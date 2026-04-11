# Get Training Indices for a Fold

An S3 generic method to retrieve the training sample indices for a
specific fold from a cross-validation specification object.

## Usage

``` r
train_indices(obj, fold_num, ...)
```

## Arguments

- obj:

  A cross-validation specification object (e.g., inheriting from
  \`cross_validation\`).

- fold_num:

  An integer specifying the fold number for which to retrieve training
  indices.

- ...:

  Additional arguments passed to methods.

## Value

An integer vector of training indices.

## Examples

``` r
cval <- kfold_cross_validation(len = 20, nfolds = 4)
train_indices(cval, 1)
#>  [1]  2  3  4  5  6  7  8  9 10 12 15 16 18 19 20
```
