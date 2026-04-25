# Constructing Datasets for MVPA Analysis

This vignette shows how to construct neuroimaging datasets for MVPA. We
start with a runnable synthetic example, then show the pattern for real
volumetric and surface data.

## Quick Start with Synthetic Data

The easiest way to get started is with
[`gen_sample_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md),
which creates a complete dataset and design in one call:

``` r
ds <- gen_sample_dataset(D = c(6, 6, 6), nobs = 80, blocks = 4, nlevels = 2)
print(ds$dataset)
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
print(ds$design)
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

The dataset bundles a 4-D `NeuroVec` (voxels × time) and a binary
`NeuroVol` mask. Two pictures make the shape concrete: a mid-axial slice
of one trial’s voxel pattern, and the trial × run grid of the design.

![Left: a mid-axial slice of one trial's voxel pattern (signed
activation, blue↔red). The synthetic mask covers the full 6x6x6 volume
here; in real data the mask would be a brain envelope. Right: each trial
coloured by run; runs split into equal-length blocks of
trials.](Constructing_Datasets_files/figure-html/synthetic-visual-1.png)

Left: a mid-axial slice of one trial’s voxel pattern (signed activation,
blue↔︎red). The synthetic mask covers the full 6x6x6 volume here; in real
data the mask would be a brain envelope. Right: each trial coloured by
run; runs split into equal-length blocks of trials.

You can also assemble the pieces manually if you want explicit control
over the design data frame:

``` r
dset <- mvpa_dataset(ds$dataset$train_data, mask = ds$dataset$mask)

design_df <- data.frame(
  Y     = ds$design$y_train,
  block = ds$design$block_var
)
mvdes <- mvpa_design(design_df, y_train = ~ Y, block_var = ~ block)
mvdes
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

## Creating a Real Volumetric (Image-Based) Dataset

The examples below use `eval = FALSE` because they reference file paths
you should replace with your own data.

This example assumes you have a 4D fMRI file (“bold.nii.gz”) and a
corresponding 3D brain mask file (“mask.nii.gz”).

``` r
library(neuroim2)

# Read the fMRI data as a NeuroVec object using neuroim2::read_vec
train_neurovec <- neuroim2::read_vec("path/to/bold.nii.gz", mode = "normal")

# Read the brain mask and create a NeuroVol object
mask_vec <- neuroim2::read_vec("path/to/mask.nii.gz", mode = "normal")
mask_vol <- NeuroVol(as.array(mask_vec), NeuroSpace(dim(mask_vec), spacing = c(2, 2, 2)))

# Create the MVPA image dataset
real_dataset <- mvpa_dataset(train_data = train_neurovec, mask = mask_vol)

# Display dataset details
print(real_dataset)
```

## Creating a Real Surface-Based Dataset

This example assumes you have cortical geometry stored in a file (e.g.,
“subject.lh.smoothwm.asc”) and a signal matrix in CSV format
(“surface_data.csv”). The signal matrix should have dimensions
corresponding to the number of vertices and the number of observations.

``` r
## remotes::install_github("bbuchsbaum/neurosurf")
library(neurosurf)

# Load the cortical geometry
geom <- read_surf_geometry("path/to/subject.lh.smoothwm.asc")

# Read the surface data from a CSV file
# The CSV should not have a header and have dimensions: number of vertices x number of observations
data_matrix <- as.matrix(read.csv("path/to/surface_data.csv", header = FALSE))

# Verify that the number of rows in the data matches the geometry
nvert <- nrow(neurosurf::vertices(geom))
if(nrow(data_matrix) != nvert) {
  stop("The number of vertices in the data does not match the geometry.")
}

# Create a NeuroSurfaceVector using the geometry and the data matrix
real_neurosurf <- NeuroSurfaceVector(geom, 1:nvert, data_matrix)

# Create the MVPA surface dataset; if no mask is provided, one is generated automatically
real_surface_dataset <- mvpa_surface_dataset(train_data = real_neurosurf, name = "lh")

# Display dataset details
print(real_surface_dataset)
```

## Notes

- For volumetric datasets, the NeuroSpace object is used to define
  dimensions and voxel spacing.
- For surface datasets, ensure that the signal matrix and cortical
  geometry are compatible.

Replace the file paths with actual paths to your real data files before
running the code.
