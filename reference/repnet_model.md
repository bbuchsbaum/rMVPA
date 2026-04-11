# Representational connectivity model (ReNA-RC)

For each ROI or searchlight, computes representational connectivity
between the ROI RDM and a seed RDM, optionally controlling for confound
RDMs (block, lag, behavior, etc.).

## Usage

``` r
repnet_model(
  dataset,
  design,
  repnet_des,
  distfun = cordist(method = "pearson"),
  simfun = c("pearson", "spearman"),
  ...
)
```

## Arguments

- dataset:

  An `mvpa_dataset`.

- design:

  An `mvpa_design`.

- repnet_des:

  Output of
  [`repnet_design`](http://bbuchsbaum.github.io/rMVPA/reference/repnet_design.md).

- distfun:

  Distance function for ROI RDM (e.g. `cordist(method = "pearson")`).

- simfun:

  Character: similarity metric between ROI and seed RDMs (`"pearson"` or
  `"spearman"`).

- ...:

  Extra fields stored on the model spec.

## Value

A model spec of class `"repnet_model"` compatible with
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
and
[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires repnet_design with seed_rdm
  # model <- repnet_model(dataset, design, repnet_des)
} # }
```
