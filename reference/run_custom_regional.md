# Run a Custom Analysis Function Regionally

Applies a user-defined function to the data within each specified region
of interest (ROI) and returns the results as a tibble.

## Usage

``` r
run_custom_regional(
  dataset,
  region_mask,
  custom_func,
  ...,
  .cores = 1,
  .verbose = FALSE
)
```

## Arguments

- dataset:

  An \`mvpa_dataset\` or \`mvpa_surface_dataset\` object.

- region_mask:

  A \`NeuroVol\` or \`NeuroSurface\` object where each region is
  identified by a unique integer greater than 0.

- custom_func:

  A function to apply to each ROI's data. It should accept two
  arguments:

  - \`roi_data\`: A matrix or tibble containing the data (samples x
    features) for the current ROI.

  - \`roi_info\`: A list containing \`id\` (the region number) and
    \`indices\` (the feature indices for this ROI).

  The function \*must\* return a named list or a single-row data frame
  (or tibble) containing scalar metric values.

- ...:

  Optional arguments passed to \`mvpa_iterate\` (e.g., \`batch_size\`).

- .cores:

  Number of cores to use for parallel processing via the \`future\`
  framework. Defaults to 1 (sequential). Set using \`future::plan()\`
  beforehand for more control.

- .verbose:

  Logical. If \`TRUE\`, prints progress messages during iteration.
  Defaults to \`FALSE\`.

## Value

A \`tibble\` where each row corresponds to an ROI. It includes:

- \`id\`: The ROI identifier (region number).

- Columns corresponding to the names returned by \`custom_func\`.

- \`error\`: Logical indicating if an error occurred for this ROI.

- \`error_message\`: The error message if an error occurred.

## Details

This function provides a simplified interface for applying custom
analyses per ROI without needing to define a full \`mvpa_model\`
specification or implement S3 methods. It leverages the parallel
processing and iteration capabilities of \`rMVPA\`.

The user-supplied \`custom_func\` performs the core calculation for each
ROI. The framework handles extracting data, iterating over ROIs
(potentially in parallel), catching errors from \`custom_func\`, and
formatting the output into a convenient flat table.

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate sample dataset
dset_info <- gen_sample_dataset(D = c(8,8,8), nobs = 50, nlevels = 2)
dataset_obj <- dset_info$dataset

# Create a region mask with 3 ROIs
mask_arr <- array(0, dim(dataset_obj$mask))
mask_arr[1:4, 1:4, 1:4] <- 1
mask_arr[5:8, 1:4, 1:4] <- 2
mask_arr[1:4, 5:8, 5:8] <- 3
region_mask_vol <- neuroim2::NeuroVol(mask_arr, neuroim2::space(dataset_obj$mask))

# Define a custom function: calculate mean and sd for each ROI
my_roi_stats <- function(roi_data, roi_info) {
  mean_signal <- mean(roi_data, na.rm = TRUE)
  sd_signal <- sd(roi_data, na.rm = TRUE)
  list(mean_signal = mean_signal, sd_signal = sd_signal, n_features = ncol(roi_data))
}

custom_results <- run_custom_regional(dataset_obj, region_mask_vol, my_roi_stats)
print(custom_results)
} # }
```
