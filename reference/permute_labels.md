# Permute Training Labels in an MVPA Design

Returns a copy of `design` with `y_train`, `cv_labels`, and `targets`
replaced by permuted versions. `block_var` is never permuted.

## Usage

``` r
permute_labels(
  design,
  method = c("within_block", "circular_shift", "global"),
  seed = NULL
)
```

## Arguments

- design:

  An `mvpa_design` object.

- method:

  Character. One of `"within_block"`, `"circular_shift"`, or `"global"`.

- seed:

  Optional integer seed; RNG state is restored on exit.

## Value

A modified `mvpa_design` with permuted labels.
