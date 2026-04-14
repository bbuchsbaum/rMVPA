# Process ROI

Process a region of interest (ROI) and return the formatted results.

## Usage

``` r
process_roi(mod_spec, roi, rnum, ...)

# Default S3 method
process_roi(mod_spec, roi, rnum, center_global_id = NA, ...)

# S3 method for class 'custom_internal_model_spec'
process_roi(mod_spec, roi, rnum, ...)
```

## Arguments

- mod_spec:

  The model specification object.

- roi:

  The region of interest data.

- rnum:

  A numeric or string identifier for the ROI.

- ...:

  Additional arguments passed to the method-specific function.

- center_global_id:

  Optional global ID of the center voxel. Defaults to NA.

## Value

A tibble row containing the performance metrics for the ROI.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(4, 4, 4), 20, blocks = 2)
  cv <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  spec <- mvpa_model(
    model = mdl,
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
  vox <- sample(which(ds$dataset$mask > 0), 30)
  samp <- data_sample(ds$dataset, vox)
  roi_obj <- as_roi(samp, ds$dataset)
  process_roi(spec, roi_obj, 1)
#> ERROR [2026-04-14 02:18:01] Class model prediction failed: length of 'dimnames' [2] not equal to array extent
#> WARN [2026-04-14 02:18:01] error fitting model 1 : Prediction failed: length of 'dimnames' [2] not equal to array extent
#> # A tibble: 1 × 8
#>   result indices performance    id error error_message   warning warning_message
#>   <list> <list>  <list>      <dbl> <lgl> <chr>           <lgl>   <chr>          
#> 1 <NULL> <dbl>   <NULL>          1 TRUE  Prediction fai… FALSE   ~              
# }
```
