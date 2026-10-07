# MS-ReVE: build temporal nuisance RDMs from a spec

Convenience helper that produces a named list of condition-level
temporal nuisance RDMs for MS-ReVE analyses given an `mvpa_design`, a
time index (trial onsets or order), and a compact specification. Each
entry can be a kernel-based nuisance (delegates to
`temporal_nuisance_for_msreve`) or an HRF-overlap-based nuisance built
trial-wise then reduced to condition level.

## Usage

``` r
msreve_temporal_confounds(
  mvpa_design,
  time_idx,
  spec,
  reduce = c("min", "mean", "median", "nn_min"),
  within_blocks_only = TRUE
)
```

## Arguments

- mvpa_design:

  an `mvpa_design` object

- time_idx:

  numeric/integer time index per trial (seconds, TRs, or ordinal)

- spec:

  a named list of entries. For kernel entries: include `kernel` and
  optional `width,power,lambda,sigma,units,TR,metric`. For HRF entries,
  set `kind="hrf"` and include `TR` plus optional
  `durations,hrf,oversampling,length_s,similarity,metric`.

- reduce:

  reduction method for condition-level aggregation
  (`"min"`,`"mean"`,`"median"`,`"nn_min"`)

- within_blocks_only:

  logical; if TRUE ignore across-run pairs

## Value

a named list of K x K matrices aligned to levels(mvpa_design\$Y)

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 4, blocks = 3)
spec <- list(lag = list(kernel = "exp", lambda = 3),
             hrf = list(kind = "hrf", TR = 2))
conf <- msreve_temporal_confounds(ds$design,
                                  time_idx = seq(0, by = 2, length.out = 24),
                                  spec = spec)
lapply(conf, dim)
#> $lag
#> [1] 4 4
#> 
#> $hrf
#> [1] 4 4
#> 
```
