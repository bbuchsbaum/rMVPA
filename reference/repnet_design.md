# Representational connectivity design helper (ReNA-RC)

Build a design object for representational connectivity, specifying item
keys, a seed RDM, and optional confound RDMs.

## Usage

``` r
repnet_design(design, key_var, seed_rdm, confound_rdms = NULL)
```

## Arguments

- design:

  An `mvpa_design` object (as used elsewhere in rMVPA).

- key_var:

  Column name or formula giving item identity (e.g. `~ ImageID`).

- seed_rdm:

  A K x K matrix or `"dist"` object; rows/cols labelled by item IDs.

- confound_rdms:

  Optional named list of K x K matrices or `"dist"` objects, used as
  nuisance RDMs (e.g. block, lag, behavior).

## Value

A list with fields:

- `key`: factor of item IDs (length = nrow(design\$train_design))

- `seed_rdm`: seed RDM as a matrix with row/colnames

- `confound_rdms`: named list of confound RDM matrices

## Examples

``` r
trials <- data.frame(ImageID = rep(letters[1:4], 3), run = rep(1:3, each = 4))
design <- mvpa_design(trials, y_train = ~ ImageID, block_var = ~ run)
my_rdm <- as.matrix(dist(matrix(rnorm(4 * 3), 4, 3)))
rownames(my_rdm) <- colnames(my_rdm) <- letters[1:4]
des <- repnet_design(design, ~ ImageID, seed_rdm = my_rdm)
str(des, max.level = 1)
#> List of 3
#>  $ key          : Factor w/ 4 levels "a","b","c","d": 1 2 3 4 1 2 3 4 1 2 ...
#>  $ seed_rdm     : num [1:4, 1:4] 0 3.04 2.23 1.79 3.04 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ confound_rdms: list()
```
