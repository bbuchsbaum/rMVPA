# Get Started with rMVPA

## Why rMVPA?

Decoding stimulus categories from fMRI activation patterns requires
coordinating data loading, model selection, cross-validation, and
spatial iteration (searchlights or ROIs). rMVPA handles all of this so
you can focus on your scientific question rather than pipeline plumbing.

This vignette walks you through three analyses on synthetic data:

1.  **Searchlight classification** – decode conditions at every voxel
    neighbourhood in the brain.
2.  **Regional classification** – decode within predefined ROIs.
3.  **RSA (Representational Similarity Analysis)** – test whether neural
    patterns match a predicted similarity structure.

## Quick example: searchlight classification

A complete searchlight run is four short steps. Each builds one object
the next step needs.

**1. Load the package and generate a tiny synthetic dataset.**
[`gen_sample_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md)
returns both the imaging data (an `mvpa_dataset`) and the experimental
design (an `mvpa_design`).

``` r

library(rMVPA)
library(neuroim2)

data <- gen_sample_dataset(D = c(6, 6, 6), nobs = 80, blocks = 4,
                           nlevels = 2, response_type = "categorical",
                           data_mode = "image")
```

**2. Wrap the data in an `mvpa_dataset` and choose a cross-validation
scheme.** Run-blocked CV uses each scanning run as a held-out fold — the
standard way to keep training and test trials temporally independent.

``` r

dset <- mvpa_dataset(data$dataset$train_data, mask = data$dataset$mask)
cval <- blocked_cross_validation(data$design$block_var)
```

**3. Pick a classifier and bind it into a model specification.**
`load_model("sda_notune")` retrieves a pre-registered shrinkage
discriminant model from the `MVPAModels` registry;
[`mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md)
glues the dataset, design, classifier, and CV scheme together.

``` r

mod  <- load_model("sda_notune")
mspec <- mvpa_model(mod, dataset = dset, design = data$design,
                    crossval = cval,
                    tune_grid = data.frame(lambda = 0.01, diagonal = FALSE))
```

**4. Run a 4 mm searchlight.** Each voxel becomes the centre of a
sphere; the model is fit and cross-validated locally; performance
metrics are returned per centre as voxel-wise maps.

``` r

sl_result <- run_searchlight(mspec, radius = 4, method = "standard")
names(sl_result$results)
#> [1] "Accuracy" "AUC"
```

[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
returned a named list of performance maps (one per metric — accuracy and
AUC by default for two-class problems). Each entry is a `NeuroVol` you
can save with
[`neuroim2::write_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/write_vol-methods.html)
or display with any volumetric viewer.

## Regional classification

The same model spec runs over predefined ROIs instead of every voxel.
The only new piece is a *region mask*: an integer-labelled `NeuroVol`
where each non-zero value identifies one ROI.

**1. Build a toy 3-region mask from the active voxels.** Real analyses
use an atlas or a parcellation; here we just split the mask randomly so
the example is self-contained.

``` r

mask <- data$dataset$mask
set.seed(42)
region_mask <- NeuroVol(
  sample(1:3, sum(mask), replace = TRUE),
  space(mask),
  indices = which(mask > 0)
)
```

**2. Run regional MVPA.**
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
reuses the model spec from the searchlight example — same classifier,
same CV scheme, just a different spatial unit.

``` r

reg_result <- run_regional(mspec, region_mask)
reg_result$performance_table
#> # A tibble: 3 × 3
#>   roinum Accuracy     AUC
#>    <int>    <dbl>   <dbl>
#> 1      1    0.512 -0.0387
#> 2      2    0.512  0.151 
#> 3      3    0.462 -0.166
```

Each row reports one ROI’s cross-validated metrics. A note on the
values: rMVPA reports **AUC as `AUC − 0.5`**, so the chance level is
**0**, not 0.5, and slightly negative AUCs (e.g. `-0.04`) just mean the
classifier did marginally worse than chance — typical for a synthetic
null example like this one. Real datasets with signal will give clearly
positive AUCs.

## RSA in 30 seconds

RSA asks whether the pattern of neural similarities across conditions
matches a model-predicted similarity structure.

``` r

# Create a synthetic dataset with 5 conditions
rsa_data <- gen_sample_dataset(D = c(6, 6, 6), nobs = 100, blocks = 5,
                               nlevels = 5, response_type = "categorical",
                               data_mode = "image")

# Hypothetical model: conditions are ordered, so nearby conditions are similar
ncond <- 5
model_rdm <- as.matrix(dist(1:ncond))

# Build the RSA design and model
# data= is a list of predictor RDMs; block_var is the run vector
rsa_des <- rsa_design(~ model_rdm,
                      data = list(model_rdm = model_rdm),
                      block_var = rsa_data$design$block_var)

dset_rsa <- mvpa_dataset(rsa_data$dataset$train_data, mask = rsa_data$dataset$mask)
rsa_mod  <- rsa_model(dataset = dset_rsa, design = rsa_des,
                      distmethod = "spearman", regtype = "pearson")

# Regional RSA
rsa_result <- run_regional(rsa_mod, region_mask = rsa_data$dataset$mask)
head(rsa_result$performance_table)
#> # A tibble: 1 × 2
#>   roinum model_rdm
#>    <int>     <dbl>
#> 1      1    0.0121
```

## Where to go next

Pick the next vignette by what you want to do:

**Build a real classification analysis**

- Cross-validation that respects fMRI run structure:
  [`vignette("CrossValidation")`](http://bbuchsbaum.github.io/rMVPA/articles/CrossValidation.md)
- Searchlight pipelines, randomized variants, and result handling:
  [`vignette("Searchlight_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Searchlight_Analysis.md)
- ROI / parcellation analyses with pooled diagnostics:
  [`vignette("Regional_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Regional_Analysis.md)
- Feature selection inside the CV loop (avoid leakage):
  [`vignette("FeatureSelection")`](http://bbuchsbaum.github.io/rMVPA/articles/FeatureSelection.md)

**Test representational hypotheses**

- Standard RSA:
  [`vignette("RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/RSA.md)
- Decompose geometry by signed contrasts (MS-ReVE):
  [`vignette("Contrast_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Contrast_RSA.md)
- ROI-to-ROI connectivity through model RDMs:
  [`vignette("Model_Space_Connectivity")`](http://bbuchsbaum.github.io/rMVPA/articles/Model_Space_Connectivity.md)
- Predict neural patterns from a feature matrix:
  [`vignette("Feature_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA.md)
  (with
  [`vignette("Feature_RSA_Advanced_Workflows")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Advanced_Workflows.md)
  for extensions)
- Per-trial RSA with built-in across-block masking:
  [`vignette("Vector_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Vector_RSA.md)

**Cross-domain transfer (encoding ↔︎ retrieval, perception ↔︎ memory)**

- The naive baseline first:
  [`vignette("Naive_Cross_Decoding")`](http://bbuchsbaum.github.io/rMVPA/articles/Naive_Cross_Decoding.md)
  — also defines the ReNA / REMAP terminology used by the rest of that
  section.

**Run from the command line, not R**

- Installable wrappers and shared workflow flags:
  [`vignette("CommandLine")`](http://bbuchsbaum.github.io/rMVPA/articles/CommandLine.md)

**Plug in your own analysis**

- Wrap any per-ROI / per-sphere function:
  [`vignette("CustomAnalyses")`](http://bbuchsbaum.github.io/rMVPA/articles/CustomAnalyses.md)
- Build a full S3 plugin with metric schemas:
  [`vignette("Plugin_Development")`](http://bbuchsbaum.github.io/rMVPA/articles/Plugin_Development.md)
