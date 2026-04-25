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

``` r
library(rMVPA)
library(neuroim2)

# Generate a synthetic 6 x 6 x 6 dataset with 80 observations in 4 runs
data <- gen_sample_dataset(D = c(6, 6, 6), nobs = 80, blocks = 4,
                           nlevels = 2, response_type = "categorical",
                           data_mode = "image")

# Build the MVPA model: dataset + design + classifier + cross-validation
dset <- mvpa_dataset(data$dataset$train_data, mask = data$dataset$mask)
cval <- blocked_cross_validation(data$design$block_var)
mod  <- load_model("sda_notune")

mspec <- mvpa_model(mod, dataset = dset, design = data$design,
                    crossval = cval,
                    tune_grid = data.frame(lambda = 0.01, diagonal = FALSE))

# Run a searchlight with a 4-mm radius
sl_result <- run_searchlight(mspec, radius = 4, method = "standard")
names(sl_result$results)
#> [1] "Accuracy" "AUC"
```

[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
returns a named list of performance maps (one per metric). Each map is a
`NeuroVol` you can write to disk with
[`neuroim2::write_vol()`](https://bbuchsbaum.github.io/neuroim2/reference/write_vol-methods.html).

## Regional classification

The same model spec works for ROI-based analysis. You just need a region
mask – an integer-labeled volume where each non-zero value defines an
ROI.

``` r
# Build a toy 3-region mask from the active voxels
mask   <- data$dataset$mask
nvox   <- sum(mask)
set.seed(42)
region_mask <- NeuroVol(sample(1:3, nvox, replace = TRUE),
                        space(mask), indices = which(mask > 0))

reg_result <- run_regional(mspec, region_mask)
print(reg_result$performance_table)
#> # A tibble: 3 × 3
#>   roinum Accuracy     AUC
#>    <int>    <dbl>   <dbl>
#> 1      1    0.512 -0.0387
#> 2      2    0.512  0.151 
#> 3      3    0.462 -0.166
```

Each row of the performance table reports cross-validated metrics for
one ROI.

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
- Specialized variants and the workflow decision map:
  [`vignette("Advanced_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Advanced_RSA.md)
  and
  [`vignette("Feature_RSA_Advanced_Workflows")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Advanced_Workflows.md)

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
