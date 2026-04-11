
<!-- badges: start -->

[![R-CMD-check](https://github.com/bbuchsbaum/rMVPA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/rMVPA/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- badges: end -->

## Multivoxel Pattern Analysis in R

### This package is under development.

### Introduction

`rMVPA` is an R library for multivariate pattern analysis of
neuroimaging data. The goal of this library is to make MVPA analyses
easy. It can be used both programmatically from within R or using a
command line interface. `rMVPA` provides a lightweight model registry
and efficient resampling methods for machine learning. What `rMVPA`
provides is the infrastructure for conducting machine learning analyses
on neuroimaging data.

Documentation and vignettes: <https://bbuchsbaum.github.io/rMVPA/>

### API lifecycle

The supported public workflow/API surface is now intentionally smaller and
explicit. Use `rmvpa_api_lifecycle()` to inspect lifecycle tiers and
`rmvpa_stable_api()` to list the stable entry points intended for scripts,
extensions, and downstream packages.

### Installation

#### Using devtools

To install `rMVPA` from within R, use the `devtools` function
`install_github`. You will need the development version of `neuroim2` as
well.

From within R:

    library(devtools)
    install_github("bbuchsbaum/neuroim2")
    install_github("bbuchsbaum/rMVPA")

#### Using git from the command line

    git clone git@github.com:bbuchsbaum/rMVPA.git
    R CMD install rMVPA

### Command line scripts

`rMVPA` now ships packaged CLI wrappers with native command names:

- `rmvpa-searchlight`
- `rmvpa-regional`

Install the package:

```r
pak::pak("bbuchsbaum/rMVPA")
```

Then copy the wrappers into a directory on your `PATH`:

```r
rMVPA::install_cli("~/.local/bin", overwrite = TRUE)
```

If needed, add that directory to `PATH`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

On macOS/Linux, you can then run:

```bash
rmvpa-searchlight --help
rmvpa-regional --help
rmvpa-searchlight --version
```

The unified CLI vignette is `vignette("CommandLine")`. The legacy source-tree
entrypoints remain available as compatibility shims, but the packaged commands
above are the supported interface.

## Citation

If you use rMVPA in your research, please cite:

Buchsbaum, B. (2026). *rMVPA: Multivoxel Pattern Analysis in R*. R package version 0.1.2. <https://github.com/bbuchsbaum/rMVPA>
