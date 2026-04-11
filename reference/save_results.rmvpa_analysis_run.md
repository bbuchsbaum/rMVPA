# Save a High-Level Analysis Run

Writes all results from an `rmvpa_analysis_run` and a root manifest
describing the workflow configuration and per-entry outputs.

## Usage

``` r
# S3 method for class 'rmvpa_analysis_run'
save_results(
  x,
  dir = NULL,
  level = c("standard", "minimal", "complete"),
  stack = c("none", "auto", "vec"),
  fname = "analysis.nii.gz",
  include = NULL,
  dtype = NULL,
  overwrite = FALSE,
  quiet = FALSE
)
```

## Arguments

- x:

  An `rmvpa_analysis_run`.

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

Invisibly, a nested list of file paths.
