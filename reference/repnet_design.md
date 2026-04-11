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
if (FALSE) { # \dontrun{
  des <- repnet_design(design, ~ ImageID, seed_rdm=my_rdm)
} # }
```
