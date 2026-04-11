# Build multiple temporal confounds from a spec

Convenience function to produce a named list of temporal nuisance RDMs
from a compact specification. Items can be kernel-based (`temporal_rdm`)
or HRF-overlap based (`temporal_hrf_overlap`).

## Usage

``` r
temporal_confounds(
  spec,
  onsets,
  run = NULL,
  units = c("auto", "sec", "TR", "index"),
  TR = NULL,
  as_dist = TRUE
)
```

## Arguments

- spec:

  a named list. Each element is a list describing one confound. For
  kernel-based items, include at least `kernel`. For HRF items, set
  `kind="hrf"` and provide `TR`. Common fields: `within_blocks_only`,
  `normalize`, `metric`.

- onsets:

  numeric vector of event onsets (seconds or TRs depending on `units`)

- run:

  optional vector of run/block identifiers

- units:

  one of `"auto"`, `"sec"`, `"TR"`, or `"index"`

- TR:

  repetition time in seconds (optional for kernel, required for HRF)

- as_dist:

  logical; if TRUE return `dist` objects (default TRUE)

## Value

A named list of `dist` objects (or matrices if `as_dist=FALSE`)

## Examples

``` r
if (FALSE) { # \dontrun{
  spec <- list(lag=list(kernel="exp", lambda=3), hrf=list(kind="hrf"))
  conf <- temporal_confounds(spec, onsets=1:20, run=rep(1:4,each=5), TR=2)
} # }
```
