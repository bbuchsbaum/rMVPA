# Temporal RDM wrapper for formula usage

Convenience wrapper for
[`temporal_rdm`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md)
that simplifies usage in RSA formulas.

## Usage

``` r
temporal(index, block = NULL, ..., as_dist = TRUE)
```

## Arguments

- index:

  numeric or integer vector representing trial order or time

- block:

  optional vector of run/block identifiers

- ...:

  additional parameters passed to
  [`temporal_rdm`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md)

- as_dist:

  logical; if TRUE return a dist object (default TRUE)

## Value

A dist object or matrix representing temporal relationships

## Details

This function provides a shorter name for use in RSA design formulas. It
calls `temporal_rdm` with the same parameters.

## See also

[`temporal_rdm`](https://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md)

## Examples

``` r
run_ids <- rep(1:4, each = 10)
set.seed(1)
task_rdm <- dist(matrix(rnorm(40 * 3), 40, 3))

# Build a temporal nuisance RDM and use it in an RSA design
rdes <- rsa_design(
  ~ task_rdm + temp_rdm,
  data = list(task_rdm = task_rdm,
              temp_rdm = temporal(1:40, block = run_ids, kernel = "adjacent", width = 2),
              run = run_ids),
  block_var = ~ run
)
```
