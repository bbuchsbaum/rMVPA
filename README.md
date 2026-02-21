
<!-- badges: start -->

[![R-CMD-check](https://github.com/bbuchsbaum/rMVPA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/rMVPA/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

## Multivoxel Pattern Analysis in R

### This package is under development.

### Introduction

`rMVPA` is an R library for multivariate pattern analysis of
neuroimaging data. The goal of this library is to make MVPA analyses
easy. It can be used both programmatically from within R or using a
command line interface. ‘rMVPA’ provides a lightweight model registry
and efficient resampling methods for machine learning. What `rMVPA`
provides is the infrastructure for conducting machine learning analyses
on neuroimaging data.

Documentation and vignettes: <https://bbuchsbaum.github.io/rMVPA/>

Method note: ITEM vs continuous-time hrfdecoder paths:
`inst/notes/item_vs_hrfdecoder.md`

### Searchlight Performance Defaults

Searchlight execution now uses a deterministic optimized runtime policy by
default (no ambient option/profile/env switches):

- fold-cache reuse
- geometry cache reuse
- clustered-neighbor fastpath
- matrix-first ROI processing
- fast ROI filtering
- RSA fast kernel
- naive cross-decoding fast kernel
- backend policy default = `auto`
- SWIFT multiclass searchlight engine (eligible multiclass, standard CV paths; outputs `Accuracy`, `AUC`, and `SWIFT_Info`)

To force a specific backend for a call, pass `backend=` explicitly:

```r
run_searchlight(mspec, radius = 8, method = "standard", backend = "default")
```

To force a specific searchlight engine for a call, pass `engine=` explicitly:

```r
run_searchlight(mspec, radius = 8, method = "standard", engine = "legacy")
```

Engine selection is recorded on results for auditability:

- `attr(res, "searchlight_engine")`
- `attr(res, "searchlight_engine_requested")`
- `attr(res, "searchlight_engine_resolved")`
- `attr(res, "searchlight_engine_reason")`
- `attr(res, "searchlight_engine_fallback_error")` (only when fallback occurred)

Operational rollback and release checklist:
`docs/searchlight_performance_rollbacks.md`

Perf guardrail harness (`scripts/run_perf_guardrails.R`) now enforces a
singleton lock by default to prevent concurrent local runs from polluting
timings. Override controls:

- `RMVPA_PERF_SINGLETON=true|false` (default `true`)
- `RMVPA_PERF_LOCK_DIR=/path/to/lock` (default `./.tmp/perf_guardrails.lock`)
- `RMVPA_PERF_LOCK_STALE_SECONDS=<seconds>` (default `21600`)

### Installation

### Using devtools

To install `rMVPA` from within R, use the `devtools` function
`install_github`. You will need the development version of `neuroim2` as
well.

From within R:

    #library(devtools)
    install_github("bbuchsbaum/neuroim2")
    install_github("bbuchsbaum/rMVPA")

### Using `git` from the command line

    git clone git@github.com:bbuchsbaum/rMVPA.git
    R CMD install rMVPA

### Optionally install command line scripts for “coding-free” MVPA analysis:

`wget https://raw.githubusercontent.com/bbuchsbaum/rMVPA/master/scripts/MVPA_Searchlight.R`

`wget https://raw.githubusercontent.com/bbuchsbaum/rMVPA/master/scripts/MVPA_Regional.R`

Then, move these files to a folder on your `PATH` and make them
executable:

`chmod +x MVPA_Searchlight.R`

`chmod +x MVPA_Regional.R`


## Albers theme
This package uses the albersdown theme. Vignettes are styled with `vignettes/albers.css` and a local `vignettes/albers.js`; the palette family is provided via `params$family` (default 'red'). The pkgdown site uses `template: { package: albersdown }`.


## Albers theme
This package uses the albersdown theme. Vignettes are styled with `vignettes/albers.css` and a local `vignettes/albers.js`; the palette family is provided via `params$family` (default 'red'). The pkgdown site uses `template: { package: albersdown }`.

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The palette family is provided via `params$family` (default 'red'). The pkgdown site uses `template: { package: albersdown }`.
<!-- albersdown:theme-note:end -->
