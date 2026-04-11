# Temporal nuisance RDM at condition level (for MS-ReVE)

Creates a condition-level temporal nuisance RDM for use with
MS-ReVE/contrast_rsa_model. This function reduces trial-level temporal
relationships to condition-level relationships.

## Usage

``` r
temporal_nuisance_for_msreve(
  mvpa_design,
  time_idx,
  reduce = c("min", "mean", "median", "nn_min"),
  kernel = c("adjacent", "boxcar", "linear", "poly", "exp", "gauss"),
  units = c("auto", "sec", "TR", "index"),
  TR = NULL,
  metric = c("distance", "similarity"),
  within_blocks_only = TRUE,
  ...
)
```

## Arguments

- mvpa_design:

  an mvpa_design object containing Y (condition labels) and block_var

- time_idx:

  numeric/integer vector of temporal indices, length =
  nrow(mvpa_design\$train_design)

- reduce:

  character string specifying reduction method:

  min

  :   Minimum lag between any pair of trials from two conditions

  mean

  :   Average lag between all pairs of trials

  median

  :   Median lag between all pairs of trials

  nn_min

  :   Minimum nearest-neighbor distance

- kernel:

  character string specifying kernel type (see
  [`temporal_rdm`](http://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md))

- units:

  one of `"auto"`, `"sec"`, `"TR"`, or `"index"`; used to select
  sensible defaults for decay

- TR:

  repetition time in seconds (optional)

- metric:

  return type: `"distance"` or `"similarity"` (default `"distance"`)

- within_blocks_only:

  logical; if TRUE, ignore pairs spanning different blocks (default
  TRUE)

- ...:

  additional kernel parameters passed to kernel functions (width, power,
  lambda, sigma)

## Value

K x K symmetric matrix (0 diagonal) aligned to levels(mvpa_design\$Y)

## Details

This function is designed for MS-ReVE analyses where temporal confounds
need to be modeled at the condition level rather than trial level. It
computes aggregate temporal relationships between conditions based on
the temporal structure of individual trials.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create temporal nuisance for MS-ReVE
temp_K <- temporal_nuisance_for_msreve(
  mvpa_design = mvpa_des,
  time_idx = seq_len(nrow(mvpa_des$train_design)),
  reduce = "min",
  kernel = "exp", 
  lambda = 3,
  within_blocks_only = TRUE
)

# Use in msreve_design
msreve_des <- msreve_design(
  mvpa_design = mvpa_des,
  contrast_matrix = C_mat,
  nuisance_rdms = list(temp_decay = temp_K)
)
} # }
```
