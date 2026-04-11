# Advanced Regional MVPA Analysis in rMVPA

## Introduction

Regional MVPA evaluates prediction performance within predefined brain
regions. This vignette walks through the complete workflow using rMVPA:

We’ll generate synthetic data, define ROIs with a region mask, build the
MVPA model with cross-validation, run the analysis, and examine the
results. The implementation follows the approach in
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md),
[`mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md),
and
[`mvpa_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md).

## Data Generation and Preparation

We begin by generating a synthetic volumetric dataset using the
[`gen_sample_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md)
function. This function creates a 4D array (with spatial dimensions and
multiple observations), along with a binary mask and an associated
design for cross-validation.

``` r
library(rMVPA)
library(neuroim2)
# Generate a synthetic dataset with dimensions 6x6x6, 80 observations, divided into 4 blocks
data_out <- rMVPA::gen_sample_dataset(D = c(6,6,6), nobs = 80, blocks = 4, nlevels = 2)
print(data_out)
#> $dataset
#> 
#>  MVPA Dataset 
#> 
#> - Training Data 
#>   - Dimensions:  6 x 6 x 6 x 80 observations 
#>   - Type:  DenseNeuroVec 
#> - Test Data 
#>   -  None 
#> - Mask Information 
#>   - Areas:  TRUE : 216 
#>   - Active voxels/vertices:  216 
#> 
#> 
#> $design
#> 
#>  MVPA Design 
#> 
#> - Training Data 
#>   - Observations:  80 
#>   - Response Type:  Factor
#>   - Levels:  a, b 
#>   - Class Distribution:  a: 40, b: 40 
#> - Test Data 
#>   -  None 
#> - Structure 
#>   - Blocking:  Present
#>   - Number of Blocks:  4 
#>   - Mean Block Size:  20  (SD:  0 ) 
#>   - Split Groups:  None
```

The returned list contains:

- **dataset**: an MVPA dataset object with training data and a binary
  mask.
- **design**: an MVPA design object specifying the response variable and
  block structure.

### Creating a Region Mask

For regional analysis, we need to define ROIs. Here, we create a region
mask by randomly assigning each active voxel in the binary mask a region
label (from 1 to 3). This simulates a scenario where the brain is
partitioned into three regions of interest.

``` r
# Extract the binary mask from the dataset
mask <- data_out$dataset$mask
nvox <- sum(mask)

# Create a regional mask: assign each voxel a random region number (1 to 3)
set.seed(123)  # for reproducibility
region_mask <- neuroim2::NeuroVol(sample(1:3, size = nvox, replace = TRUE), neuroim2::space(mask), indices = which(mask > 0))
table(region_mask)
#> region_mask
#>  1  2  3 
#> 69 74 73
```

## Setting Up the MVPA Model

Next, we create an MVPA model to evaluate a classification task. In
brief, we construct an `mvpa_dataset`, specify the design (including the
block variable and response), and define the model with
[`mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md)
using a chosen classifier and cross‑validation strategy.

``` r
# Create MVPA dataset object from the generated training data and mask
 dset <- mvpa_dataset(data_out$dataset$train_data, mask = data_out$dataset$mask)

# Build cross-validation structure using block information from the design
cval <- blocked_cross_validation(data_out$design$block_var)

# Load a classification model; here we use "sda" (Shrinkage Discriminant Analysis)
mod <- load_model("sda")
tune_grid <- data.frame(lambda = 0.01, diagonal = FALSE)

# Create the MVPA model object
mvpa_mod <- mvpa_model(mod, dataset = dset, design = data_out$design, crossval = cval, tune_grid = tune_grid)
print(mvpa_mod)
#> mvpa_model object. 
#> model:  sda 
#> model type:  classification 
#> tune_reps:  15 
#> tune_grid:  
#>   lambda diagonal
#> 1   0.01    FALSE
#> 
#>  Blocked Cross-Validation 
#> 
#> - Dataset Information 
#>   - Observations:  80 
#>   - Number of Folds:  4 
#> - Block Information 
#>   - Total Blocks:  4 
#>   - Mean Block Size:  20  (SD:  0 ) 
#>   - Block Sizes:  1: 20, 2: 20, 3: 20, 4: 20 
#> 
#> 
#>  MVPA Dataset 
#> 
#> - Training Data 
#>   - Dimensions:  6 x 6 x 6 x 80 observations 
#>   - Type:  DenseNeuroVec 
#> - Test Data 
#>   -  None 
#> - Mask Information 
#>   - Areas:  TRUE : 216 
#>   - Active voxels/vertices:  216 
#> 
#> 
#>  MVPA Design 
#> 
#> - Training Data 
#>   - Observations:  80 
#>   - Response Type:  Factor
#>   - Levels:  a, b 
#>   - Class Distribution:  a: 40, b: 40 
#> - Test Data 
#>   -  None 
#> - Structure 
#>   - Blocking:  Present
#>   - Number of Blocks:  4 
#>   - Mean Block Size:  20  (SD:  0 ) 
#>   - Split Groups:  None
```

The
[`mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md)
function, as defined in `mvpa_model.R`, packages all necessary
parameters including cross-validation and performance computation.

## Running the Regional Analysis

The regional analysis is executed with
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md).
Internally it prepares ROI indices from the region mask, applies the
MVPA model to each ROI, and then compiles performance metrics and
prediction tables.

``` r
# Run the regional analysis on the defined region mask
regional_results <- run_regional(mvpa_mod, region_mask)
```

The output is a `regional_mvpa_result` object with a `performance_table`
(cross‑validated metrics per region), a `prediction_table` (trial‑level
predictions), and `vol_results` (volumetric maps of performance across
the brain).

## Examining the Results

We can inspect the performance table to evaluate model accuracy in each
region.

``` r
# Display performance metrics for each region
print(regional_results$performance_table)
#> # A tibble: 3 × 3
#>   roinum Accuracy      AUC
#>    <int>    <dbl>    <dbl>
#> 1      1    0.575  0.188  
#> 2      2    0.488  0.00125
#> 3      3    0.488 -0.157
```

For a more detailed view, the prediction table shows trial-by-trial
predictions:

``` r
# Display the prediction table
print(regional_results$prediction_table)
#> # A tibble: 240 × 8
#> # Rowwise: 
#>    .rownum roinum observed pobserved predicted correct prob_a  prob_b
#>      <int>  <int> <fct>        <dbl> <fct>     <lgl>    <dbl>   <dbl>
#>  1       1      1 b           0.872  b         TRUE    0.128  0.872  
#>  2       2      1 a           0.961  a         TRUE    0.961  0.0393 
#>  3       3      1 a           0.996  a         TRUE    0.996  0.00371
#>  4       4      1 b           0.983  b         TRUE    0.0166 0.983  
#>  5       5      1 b           0.0243 a         FALSE   0.976  0.0243 
#>  6       6      1 a           0.805  a         TRUE    0.805  0.195  
#>  7       7      1 b           0.799  b         TRUE    0.201  0.799  
#>  8       8      1 a           0.0667 b         FALSE   0.0667 0.933  
#>  9       9      1 a           0.999  a         TRUE    0.999  0.00124
#> 10      10      1 a           0.554  a         TRUE    0.554  0.446  
#> # ℹ 230 more rows
```

Volumetric results (`vol_results`) can be further visualized with
neuroimaging tools to determine spatial patterns of performance.

## Under the Hood: How It Works

The
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
function internally calls
[`prep_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/prep_regional.md)
to process the region mask, and then uses
[`mvpa_iterate()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
to apply the MVPA model across each ROI. Functions such as
`combine_regional_results()` and
[`combine_prediction_tables()`](http://bbuchsbaum.github.io/rMVPA/reference/combine_prediction_tables.md)
merge the individual regional outputs into a comprehensive result.

This modular design ensures that:

- Data integrity is maintained.
- Cross-validation is properly applied.
- Results are aggregated for clear interpretation at the regional level.

## Summary

This vignette showed you how to generate synthetic neuroimaging data,
define ROIs with region masks, set up MVPA models with cross-validation,
run analyses across ROIs, and interpret the performance metrics. You now
have the tools to conduct regional MVPA analyses on your own
neuroimaging data.

For further details, see
[`vignette("Searchlight_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Searchlight_Analysis.md)
for the searchlight counterpart and
[`vignette("CrossValidation")`](http://bbuchsbaum.github.io/rMVPA/articles/CrossValidation.md)
for cross-validation strategies.
