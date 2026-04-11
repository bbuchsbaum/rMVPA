# Run a single ROI through process_roi for debugging

Run a single ROI through process_roi for debugging

## Usage

``` r
run_single_roi(
  mod_spec,
  sl_win,
  processor = NULL,
  preserve_center = TRUE,
  min_voxels = 2,
  drop_probs = FALSE
)
```

## Value

A tibble row with ROI analysis results.
