# Prepare regional data for MVPA analysis

This function processes the input data and prepares the regions for MVPA
analysis by extracting voxel indices for each region of interest (ROI)
specified in the region_mask.

## Usage

``` r
prep_regional(model_spec, region_mask)
```

## Arguments

- model_spec:

  A model specification object.

- region_mask:

  A mask representing different regions in the brain image.

## Value

A list containing information about the regions for further processing:
\* allrois: A vector of unique ROI labels. \* region_vec: A vector
representation of the region_mask. \* region_set: A sorted vector of
unique ROI labels in the region_mask. \* vox_iter: A list containing
voxel indices for each ROI. \* lens: A vector containing the number of
voxels in each ROI. \* keep: A logical vector indicating if an ROI
should be kept for analysis (those with more than one voxel).

## Examples

``` r
# Create example data
sample_data <- gen_sample_dataset(c(5, 5, 5), nobs = 100, blocks = 4)

# Create a simple region mask with 3 ROIs
mask_vol <- sample_data$dataset$mask
region_mask <- neuroim2::NeuroVol(
  sample(1:3, size = sum(mask_vol > 0), replace = TRUE),
  space = neuroim2::space(mask_vol),
  indices = which(mask_vol > 0)
)

# Create a basic model spec
model_spec <- list(dataset = sample_data$dataset)

# Prepare regional data
regional_data <- prep_regional(model_spec, region_mask)
```
