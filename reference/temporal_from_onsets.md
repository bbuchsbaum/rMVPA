# Temporal RDM from onsets (sugar)

Convenience wrapper that accepts event onsets (and optional run labels)
and delegates to `temporal_rdm`. This is useful in RSA formulas.

## Usage

``` r
temporal_from_onsets(
  onsets,
  run = NULL,
  ...,
  units = c("auto", "sec", "TR", "index"),
  TR = NULL,
  as_dist = TRUE
)
```

## Arguments

- onsets:

  numeric vector of event onsets (typically seconds or TRs)

- run:

  optional vector of run/block identifiers

- ...:

  additional parameters passed to `temporal_rdm`

- units:

  units for `onsets`: one of `"auto"`, `"sec"`, `"TR"`, or `"index"`

- TR:

  repetition time in seconds (optional)

- as_dist:

  logical; if TRUE return a dist object (default TRUE)

## Value

A dist object or matrix representing temporal relationships

## Examples

``` r
if (FALSE) { # \dontrun{
  onsets <- c(0, 2.5, 5.0, 7.5, 10.0)
  temp_rdm <- temporal_from_onsets(onsets, kernel="exp", lambda=3)
} # }
```
