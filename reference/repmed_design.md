# Representational mediation design helper (ReNA-RM)

Construct a design object for representational mediation, with predictor
(X) and outcome (Y) RDMs over items and optional confound RDMs.

## Usage

``` r
repmed_design(items, X_rdm, Y_rdm, confound_rdms = NULL)
```

## Arguments

- items:

  Character vector of item IDs (keys).

- X_rdm:

  Predictor RDM (matrix or `"dist"`).

- Y_rdm:

  Outcome RDM (matrix or `"dist"`).

- confound_rdms:

  Optional named list of confound RDMs (matrix or `"dist"`).

## Value

A list with elements:

- `items`: character vector of items

- `X_rdm`, `Y_rdm`: predictor and outcome RDMs as matrices

- `confound_rdms`: named list of confound RDM matrices

## Examples

``` r
if (FALSE) { # \dontrun{
  items <- letters[1:10]
  X_rdm <- as.matrix(dist(matrix(rnorm(10*3), 10, 3)))
  Y_rdm <- as.matrix(dist(matrix(rnorm(10*3), 10, 3)))
  des <- repmed_design(items, X_rdm, Y_rdm)
} # }
```
