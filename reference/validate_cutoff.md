# Validate Feature-Selection Cutoff

Shared helper for converting `top_k`/`top_p`-style cutoffs to an integer
feature count.

## Usage

``` r
validate_cutoff(type, value, ncol)
```

## Arguments

- type:

  Cutoff type (`top_k`, `topk`, `top_p`, `topp`).

- value:

  Cutoff value.

- ncol:

  Number of features available.

## Value

Integer number of features to retain.
