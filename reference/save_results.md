# Save MVPA Results to Disk

Generic function for writing MVPA result objects (searchlight, regional,
or simple lists of maps) to a directory on disk in a reproducible
layout.

## Usage

``` r
save_results(
  x,
  dir,
  level = c("standard", "minimal", "complete"),
  stack = c("none", "auto", "vec"),
  fname = "searchlight.nii.gz",
  include = NULL,
  dtype = NULL,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- x:

  A result object (typically a `searchlight_result`,
  `regional_mvpa_result`, or a named list of `NeuroVol`/`NeuroVec` (and
  optionally `NeuroSurface`/`NeuroSurfaceVector`) objects.

- dir:

  Directory to write into. It is created if it does not exist.

- level:

  One of `"minimal"`, `"standard"`, `"complete"`:

  - `"minimal"`: write only map files (no manifest, no tables, no aux).

  - `"standard"`: maps + a manifest describing files (default).

  - `"complete"`: maps + manifest + summary tables and auxiliary objects
    when available.

- stack:

  One of `c("none","auto","vec")`.

  - `"none"`: write one NIfTI file per metric (default).

  - `"auto"`: if all volumes are compatible (same space), stack them
    along the 4th dimension; otherwise fall back to one file per metric.

  - `"vec"`: always stack into a single 4D NIfTI (error if not
    compatible).

- fname:

  Base filename when writing a 4D file (default `"searchlight.nii.gz"`).
  Only used when `stack = "vec"` or when `stack = "auto"` and stacking
  is possible.

- include:

  Character vector of extras to include; subset of
  `c("manifest","tables","aux")`. The `level` argument sets sensible
  defaults but `include` can add to these.

- dtype:

  Optional `data_type` passed to neuroim2 writers (e.g., `"FLOAT"`,
  `"DOUBLE"`).

- overwrite:

  Logical; if `FALSE` and a target file already exists, a numeric suffix
  is appended instead of overwriting.

- quiet:

  Logical; if `TRUE`, suppress progress messages.

## Value

(invisible) a list describing what was written (file paths grouped by
type such as `maps`, `surfaces`, `tables`, `aux`, and `manifest`).

## Details

`save_results()` is an S3 generic. The main methods are:

- **`save_results.searchlight_result`**: writes volumetric maps under
  `<dir>/maps` and (optionally) surface maps under `<dir>/surfaces`.
  Depending on `level`/`include`, it can also write a per-metric summary
  table (`<dir>/tables/metric_summary.csv`), auxiliary objects (e.g.
  predictions, CV folds, design) as `.rds` files under `<dir>/aux`, and
  a manifest (`manifest.yaml/json/rds`) describing all files and their
  metric names.

- **`save_results.regional_mvpa_result`**: uses the `searchlight_result`
  method to write volumetric maps in `vol_results`, then saves regional
  tables such as `performance_table` and `prediction_table` as
  tab-delimited text. For `level != "minimal"`, ROI-wise fits (if
  present) are saved under `<dir>/fits`. A manifest can be written
  summarizing ROI counts, presence of fits, and all file paths.

- **`save_results.default`**: if `x` is a named list of
  `NeuroVol`/`NeuroVec` (and/or surface) objects, it is treated as a
  lightweight searchlight-style result and handled as above. Otherwise a
  simple `result.rds` is written to `dir`.

All methods return (invisibly) a nested list of file paths that can be
used to track outputs or drive downstream packaging/publishing.

## Examples

``` r
# \donttest{
  # After running a searchlight analysis:
  #   sl_res <- run_searchlight(mspec, radius = 2)
  # Save maps and a manifest into "sl_output"
  # save_results(sl_res, "sl_output", level = "standard")

  # After running a regional analysis:
  #   reg_res <- run_regional(mspec, regionMask)
  # Save regional maps, performance table, and fits
  # save_results(reg_res, "regional_output", level = "complete")
# }
```
