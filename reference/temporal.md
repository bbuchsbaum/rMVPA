# Temporal RDM wrapper for formula usage

Convenience wrapper for
[`temporal_rdm`](http://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md)
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
  [`temporal_rdm`](http://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md)

- as_dist:

  logical; if TRUE return a dist object (default TRUE)

## Value

A dist object or matrix representing temporal relationships

## Details

This function provides a shorter name for use in RSA design formulas. It
calls `temporal_rdm` with the same parameters.

## See also

[`temporal_rdm`](http://bbuchsbaum.github.io/rMVPA/reference/temporal_rdm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Use directly in RSA formula
rdes <- rsa_design(
  ~ task_rdm + temporal(trial_index, block=run, kernel="adjacent", width=2),
  data = list(task_rdm = task_rdm, trial_index = 1:100, run = run_ids),
  block_var = ~ run
)
} # }
```
