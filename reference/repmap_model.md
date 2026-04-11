# Representational mapping model (ReNA-Map)

For each ROI/searchlight, fits a reduced-rank regression from seed
feature vectors (X) to ROI item-level patterns (Y), summarizing mapping
rank and fit.

## Usage

``` r
repmap_model(
  dataset,
  design,
  repmap_des,
  key_var,
  rank = "auto",
  max_rank = 20,
  ridge_lambda = NULL,
  ...
)
```

## Arguments

- dataset:

  An `mvpa_dataset`.

- design:

  An `mvpa_design`.

- repmap_des:

  Output of
  [`repmap_design`](http://bbuchsbaum.github.io/rMVPA/reference/repmap_design.md).

- key_var:

  Column or formula giving item identity (e.g. `~ ImageID`).

- rank:

  `"auto"` for
  [`rrpack::cv.rrr`](https://rdrr.io/pkg/rrpack/man/cv.rrr.html) rank
  selection, integer for fixed rank, or `0` for no mapping (zero map;
  not an identity transform).

- max_rank:

  Maximum rank to search.

- ridge_lambda:

  Optional ridge penalty lambda for
  [`rrpack::rrs.fit`](https://rdrr.io/pkg/rrpack/man/rrs.fit.html).

- ...:

  Extra fields stored on the model spec.

## Value

A model spec of class `"repmap_model"`.

## Details

Internally, the item-level seed features (X) and ROI patterns (Y) are
column-centered prior to reduced-rank regression. Returned voxelwise
R-squared values are in-sample and may be negative when the mapping
underperforms the mean model; this can be useful as a diagnostic rather
than an error.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires repmap_design with seed features
  # model <- repmap_model(dataset, design, repmap_des, key_var=~ImageID)
} # }
```
