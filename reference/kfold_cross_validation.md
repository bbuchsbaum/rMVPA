# kfold_cross_validation

Construct a cross-validation specification that randomly partitions the
input set into `nfolds` folds.

## Usage

``` r
kfold_cross_validation(len, nfolds = 10)
```

## Arguments

- len:

  An integer representing the number of observations.

- nfolds:

  An integer specifying the number of cross-validation folds.

## Value

An object of class "kfold_cross_validation", "cross_validation", and
"list" containing the block_var and nfolds.

## Details

This function creates a k-fold cross-validation scheme for cases where
data needs to be split into a specified number of folds for evaluation.
It returns an object of class "kfold_cross_validation",
"cross_validation", and "list".

## Examples

``` r
cval <- kfold_cross_validation(len=100, nfolds=10)
sample_data <- as.data.frame(matrix(rnorm(100*10), 100, 10))
sample_y <- rep(letters[1:5], 20)
samples <- crossval_samples(cval, data = sample_data, y = sample_y)
stopifnot(nrow(samples) == 10)
```
