# Alias for `mvpa_multibasis_dataset`

Alias for `mvpa_multibasis_dataset`

## Usage

``` r
mvpa_multibasis_image_dataset(
  train_data,
  test_data = NULL,
  mask,
  basis_count = NULL,
  ordering = c("event_major", "basis_major"),
  basis_labels = NULL
)
```

## Arguments

- train_data:

  Training data in one of the following forms:

  - A list of `NeuroVec` objects (one per basis).

  - A character vector of 4D image paths (one file per basis).

  - A single 4D `NeuroVec` or image path containing concatenated basis
    volumes, where `basis_count` specifies splitting.

- test_data:

  Optional test data in the same format as `train_data`.

- mask:

  A `NeuroVol` mask.

- basis_count:

  Number of basis functions when using a single concatenated 4D series.

- ordering:

  Ordering of volumes in concatenated series: `"event_major"` means
  `event1(b1..bk), event2(b1..bk), ...`; `"basis_major"` means
  `b1(all events), b2(all events), ...`.

- basis_labels:

  Optional character labels for basis functions.

## Value

An object of class `mvpa_multibasis_image_dataset`.

## Examples

``` r
if (FALSE) { # \dontrun{
  # See mvpa_multibasis_dataset for examples
  ds <- gen_sample_dataset(c(5,5,5), 20)
  mb <- mvpa_multibasis_image_dataset(
    list(ds$dataset$train_data, ds$dataset$train_data),
    mask = ds$dataset$mask
  )
} # }
```
