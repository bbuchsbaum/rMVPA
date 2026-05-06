# rMVPA — multivoxel pattern analysis in R

A methodologically rigorous toolbox for decoding and representational
analysis of neuroimaging data. One model spec runs *regionally* or via
*searchlight*. Cross-validation respects fMRI run structure. RSA,
MS-ReVE contrast decomposition, model-space connectivity, ERA
cross-decoding, and domain-adaptive REMAP-RRR are all first-class — no
one-off scripts in a graveyard of unmaintained lab repos.

📖 **Documentation:** <https://bbuchsbaum.github.io/rMVPA/>

## The mental model

    mvpa_dataset  ──►   mvpa_design   ──►   model_spec   ──►   engine     ──►   result
    (voxels x time      (response,         (rsa_model,        (run_regional      (regional_mvpa_result
     + mask)             block_var,         mvpa_model,        run_searchlight)   searchlight_result)
                         splits)            feature_rsa, …)

Three S3 layers — **dataset → design → model_spec** — and two engines
(`run_regional` / `run_searchlight`). The same `model_spec` works in
both. Cross-validation, parallelism (via `future`), error handling, and
per-ROI diagnostics come for free.

## What you can build with it

| You want to … | Reach for | Vignette |
|:---|:---|:---|
| Decode condition labels from local activity (8-way category, condition pair, etc.) | [`mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md) + [`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md) / [`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md) | [`Searchlight_Analysis`](https://bbuchsbaum.github.io/rMVPA/articles/Searchlight_Analysis.html), [`Regional_Analysis`](https://bbuchsbaum.github.io/rMVPA/articles/Regional_Analysis.html), [`Haxby_2001`](https://bbuchsbaum.github.io/rMVPA/articles/Haxby_2001.html) ⓡ |
| Test whether a region’s pattern geometry matches a model RDM | [`rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md) + [`rsa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md) (or [`pair_rsa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/pair_rsa_design.md)) | [`RSA`](https://bbuchsbaum.github.io/rMVPA/articles/RSA.html), [`Kriegeskorte_92_Images`](https://bbuchsbaum.github.io/rMVPA/articles/Kriegeskorte_92_Images.html) ⓡ |
| Ask which regions share representational geometry through your models | [`model_space_connectivity()`](http://bbuchsbaum.github.io/rMVPA/reference/model_space_connectivity.md) on `rsa_model(..., return_fingerprint = TRUE)` | [`Model_Space_Connectivity`](https://bbuchsbaum.github.io/rMVPA/articles/Model_Space_Connectivity.html) |
| Decompose representational structure by signed contrasts (MS-ReVE) | [`contrast_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/contrast_rsa_model.md) | [`Contrast_RSA`](https://bbuchsbaum.github.io/rMVPA/articles/Contrast_RSA.html) |
| Cross-decode between encoding & retrieval (ERA), or perception & memory | [`era_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md), [`naive_xdec_model()`](http://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md), [`remap_rrr_model()`](http://bbuchsbaum.github.io/rMVPA/reference/remap_rrr_model.md), the ReNA family | [`ERA_RSA_Cross_Decoding`](https://bbuchsbaum.github.io/rMVPA/articles/ERA_RSA_Cross_Decoding.html), [`REMAP_RRR`](https://bbuchsbaum.github.io/rMVPA/articles/REMAP_RRR.html), and the `Naive_Cross_Decoding` glossary |
| Plug in a custom per-ROI / per-sphere function | [`run_custom_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_custom_regional.md) / [`run_custom_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_custom_searchlight.md) | [`CustomAnalyses`](https://bbuchsbaum.github.io/rMVPA/articles/CustomAnalyses.html), [`Plugin_Development`](https://bbuchsbaum.github.io/rMVPA/articles/Plugin_Development.html) |

ⓡ = real public data, not synthetic.

## A 30-second taste

``` r

library(rMVPA)

# Synthetic 6 x 6 x 6, 80 trials, 4 runs, 2 categories
data <- gen_sample_dataset(D = c(6, 6, 6), nobs = 80, blocks = 4, nlevels = 2)

dset <- mvpa_dataset(data$dataset$train_data, mask = data$dataset$mask)
cval <- blocked_cross_validation(data$design$block_var)
mod  <- load_model("sda_notune")

mspec <- mvpa_model(mod, dataset = dset, design = data$design,
                    crossval = cval,
                    tune_grid = data.frame(lambda = 0.01, diagonal = FALSE))

# Per-voxel decoding map
sl <- run_searchlight(mspec, radius = 4, method = "standard")
names(sl$results)   # "Accuracy", "AUC"
```

The `Get Started` vignette walks the same example, then routes you to
classification, RSA, and cross-domain workflows.

## Real-data demonstrations

The catalogue includes two non-synthetic vignettes that reproduce
canonical findings on public data, end-to-end:

- **Kriegeskorte 2008**
  ([vignette](https://bbuchsbaum.github.io/rMVPA/articles/Kriegeskorte_92_Images.html))
  — 92 object images, 4 subjects × 2 sessions of human-IT RDMs. Recovers
  the published model-RDM ranking (monkey-IT and animacy lead; low-level
  vision lags) and shows model-space connectivity across subjects.
- **Haxby 2001**
  ([vignette](https://bbuchsbaum.github.io/rMVPA/articles/Haxby_2001.html))
  — Subject 1 ventral-temporal patterns. Reaches **91.7 %** 8-way
  category accuracy at 12.5 % chance with `sda_notune`. Reproduces the
  canonical face / house / object confusion structure.

The bundled data lives under `inst/extdata/`; the `data-raw/*.R` scripts
regenerate it from upstream sources.

## Installation

From within R:

``` r

# install.packages("pak")
pak::pak("bbuchsbaum/neuroim2")    # required, not on CRAN yet
pak::pak("bbuchsbaum/rMVPA")
```

Or with devtools:

``` r

library(devtools)
install_github("bbuchsbaum/neuroim2")
install_github("bbuchsbaum/rMVPA")
```

## Command-line interface

`rMVPA` ships packaged CLI wrappers for searchlight and regional
analyses:

``` r

# After installing, copy the wrappers onto your PATH
rMVPA::install_cli("~/.local/bin", overwrite = TRUE)
```

``` bash
rmvpa-searchlight --help
rmvpa-regional   --help
```

See
[`vignette("CommandLine")`](https://bbuchsbaum.github.io/rMVPA/articles/CommandLine.html)
for input formats, configuration files, and parallelisation.

## API stability

Public-API surface is intentionally narrow and explicit. Use
[`rmvpa_api_lifecycle()`](http://bbuchsbaum.github.io/rMVPA/reference/rmvpa_api_lifecycle.md)
to inspect lifecycle tiers and
[`rmvpa_stable_api()`](http://bbuchsbaum.github.io/rMVPA/reference/rmvpa_stable_api.md)
to list the stable entry points intended for scripts, extensions, and
downstream packages.

## Citation

If you use rMVPA in your research, please cite:

> Buchsbaum, B. (2026). *rMVPA: Multivoxel Pattern Analysis in R*. R
> package version 0.1.2. <https://github.com/bbuchsbaum/rMVPA>
