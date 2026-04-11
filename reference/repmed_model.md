# Representational mediation model (ReNA-RM)

For each ROI/searchlight, tests whether the ROI RDM mediates the
relationship between a predictor RDM (X) and an outcome RDM (Y),
optionally including confound RDMs.

## Usage

``` r
repmed_model(
  dataset,
  design,
  repmed_des,
  key_var,
  distfun = cordist(method = "pearson"),
  ...
)
```

## Arguments

- dataset:

  An `mvpa_dataset`.

- design:

  An `mvpa_design`.

- repmed_des:

  Output of
  [`repmed_design`](http://bbuchsbaum.github.io/rMVPA/reference/repmed_design.md).

- key_var:

  Column or formula giving item identity in the mvpa_design (e.g.
  `~ ImageID`).

- distfun:

  Distance function to compute the ROI RDM (e.g.
  `cordist(method = "pearson")` ).

- ...:

  Extra fields stored on the model spec.

## Value

A model spec of class `"repmed_model"`.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires repmed_design with X_rdm and Y_rdm
  # model <- repmed_model(dataset, design, repmed_des, key_var=~ImageID)
} # }
```
