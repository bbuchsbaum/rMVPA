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

## Next steps

You now have the basic building blocks. Explore the full documentation:

- [`vignette("CrossValidation")`](http://bbuchsbaum.github.io/rMVPA/articles/CrossValidation.md)
  – cross-validation strategies for fMRI
- [`vignette("Searchlight_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Searchlight_Analysis.md)
  – in-depth searchlight workflows
- [`vignette("Regional_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Regional_Analysis.md)
  – regional analysis details
- [`vignette("RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/RSA.md)
  – full RSA tutorial
- [`vignette("CustomAnalyses")`](http://bbuchsbaum.github.io/rMVPA/articles/CustomAnalyses.md)
  – plug in your own analysis functions
- [`vignette("FeatureSelection")`](http://bbuchsbaum.github.io/rMVPA/articles/FeatureSelection.md)
  – dimensionality reduction inside CV
