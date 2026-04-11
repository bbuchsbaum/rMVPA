# Generate Sample Dataset for MVPA Analysis

Creates a synthetic dataset for testing and demonstration of MVPA
analyses.

## Usage

``` r
gen_sample_dataset(
  D,
  nobs,
  response_type = c("categorical", "continuous"),
  data_mode = c("image", "surface"),
  spacing = c(1, 1, 1),
  blocks = 5,
  nlevels = 5,
  external_test = FALSE,
  ntest_obs = nobs,
  split_by = NULL,
  na_cols = 0
)
```

## Arguments

- D:

  The data dimension(s): vector of length 2 or 3 for image data, or
  single number for surface data

- nobs:

  The number of observations

- response_type:

  Either 'categorical' or 'continuous'

- data_mode:

  Either 'image' or 'surface'

- spacing:

  The voxel spacing (default: c(1,1,1))

- blocks:

  The number of 'blocks' in the data (for cross-validation)

- nlevels:

  The number of category levels (only used if
  response_type='categorical')

- external_test:

  Whether to generate an external test set

- ntest_obs:

  The number of test observations (default: nobs)

- split_by:

  Optional factor for splitting analyses

- na_cols:

  The number of columns to randomly set to NA (default: 0)

## Value

A list containing:

- dataset:

  An `mvpa_dataset` object containing:

  - `train_data`: Training data as `NeuroVec` or `ROISurface`

  - `test_data`: Test data (if external_test=TRUE)

  - `mask`: Binary mask indicating valid voxels/vertices

- design:

  An `mvpa_design` object containing:

  - `y_train`: Response variable for training

  - `y_test`: Response variable for test set (if external_test=TRUE)

  - `block_var`: Block variable for cross-validation

  - `split_by`: Optional splitting factor

## Examples

``` r
# Generate categorical image dataset
dataset <- gen_sample_dataset(
  D = c(10,10,10),
  nobs = 100,
  response_type = "categorical",
  data_mode = "image",
  blocks = 3,
  nlevels = 2
)

# Generate continuous surface dataset
surf_data <- gen_sample_dataset(
  D = 1000,  # number of vertices
  nobs = 50,
  response_type = "continuous",
  data_mode = "surface"
)
#> Warning: RGL: unable to open X11 display
#> Warning: 'rgl.init' failed, will use the null device.
#> See '?rgl.useNULL' for ways to avoid this warning.
#> loading /home/runner/work/_temp/Library/neurosurf/extdata/std.8_lh.inflated.asc

# Generate dataset with external test set
test_dataset <- gen_sample_dataset(
  D = c(8,8,8),
  nobs = 80,
  response_type = "categorical",
  nlevels = 3,
  external_test = TRUE
)
#> external test
```
