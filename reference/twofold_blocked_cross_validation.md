# twofold_blocked_cross_validation

Construct a cross-validation specification that randomly partitions the
input set into two sets of blocks.

## Usage

``` r
twofold_blocked_cross_validation(block_var, nreps = 10)
```

## Arguments

- block_var:

  An integer vector representing the cross-validation blocks. Each block
  is indicated by a unique integer.

- nreps:

  An integer specifying the number of repetitions for the twofold split.

## Value

An object of class "twofold_blocked_cross_validation",
"cross_validation", and "list" containing the block_var, nfolds (fixed
at 2 for this function), nreps, and block_ind.

## Details

This function creates a cross-validation scheme for cases where data is
organized into blocks, and these blocks are divided into two groups for
evaluation. This approach can be useful when there is an inherent
structure or dependency within the blocks, and separating them can help
to avoid biased estimates of model performance. It returns an object of
class "twofold_blocked_cross_validation", "cross_validation", and
"list".

## Examples

``` r
blockvar <- rep(1:5, each=10)
nreps <- 5
cval <- twofold_blocked_cross_validation(blockvar, nreps=nreps)
samples <- crossval_samples(cval, as.data.frame(matrix(rnorm(50*50),50,50)), y=rep(letters[1:5],10))
stopifnot(nrow(samples) == nreps)
```
