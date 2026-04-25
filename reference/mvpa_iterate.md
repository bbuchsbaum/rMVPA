# Iterate MVPA Analysis Over Multiple ROIs

Performs multivariate pattern analysis (MVPA) across multiple regions of
interest (ROIs) using batch processing and parallel computation.

## Usage

``` r
mvpa_iterate(
  mod_spec,
  vox_list,
  ids = 1:length(vox_list),
  batch_size = NULL,
  verbose = TRUE,
  processor = NULL,
  analysis_type = c("searchlight", "regional"),
  drop_probs = FALSE,
  fail_fast = FALSE,
  save_rdm_vectors_dir = NULL
)
```

## Arguments

- mod_spec:

  An MVPA model specification object containing the dataset to analyze,
  compute_performance (logical indicating whether to compute performance
  metrics), and return_predictions (logical indicating whether to return
  predictions).

- vox_list:

  A list of voxel indices or coordinates defining each ROI to analyze.

- ids:

  Vector of identifiers for each ROI analysis. Defaults to
  1:length(vox_list).

- batch_size:

  Integer specifying number of ROIs to process per batch. For
  searchlight analyses the default is 10% of total ROIs. For regional
  analyses the default is all ROIs in one batch, unless the estimated
  extraction memory would exceed the budget set by
  `options(rMVPA.regional_mem_budget)` (default 2 GB), in which case
  batches are automatically sized to stay within the budget.

- verbose:

  Logical indicating whether to print progress messages. Defaults to
  TRUE. When `TRUE` and the progressr package is installed, a real-time
  progress bar is shown that updates as each ROI completes – even when
  running on parallel future workers. Without progressr, only coarse
  batch-level log messages are printed. Install with
  `install.packages("progressr")` and activate once per session with
  `progressr::handlers(global = TRUE)`.

- processor:

  Optional custom processing function. If NULL, uses default processor.
  Must accept parameters (obj, roi, rnum) and return a tibble.

- analysis_type:

  Character indicating the type of analysis. Defaults to "searchlight".

- drop_probs:

  Logical; if TRUE, drop per-ROI probability matrices after computing
  metrics. Default FALSE.

- fail_fast:

  Logical; if TRUE, stop immediately on first ROI error. Default FALSE.

- save_rdm_vectors_dir:

  Optional directory for writing file-backed feature-RSA RDM vector
  batches instead of retaining all vectors in memory.

## Value

A tibble containing results for each ROI with columns:

- result:

  List column of analysis results (NULL if return_predictions=FALSE).

- indices:

  List column of ROI indices used.

- performance:

  List column of performance metrics (if computed).

- id:

  ROI identifier.

- error:

  Logical indicating if an error occurred.

- error_message:

  Error message if applicable.

- warning:

  Logical indicating if a warning occurred.

- warning_message:

  Warning message if applicable.

## Details

The function processes ROIs in batches to manage memory usage. For each
batch:

1.  Extracts ROI data from the dataset.

2.  Filters out ROIs with fewer than 2 voxels.

3.  Processes each ROI using either the default or custom processor.

4.  Combines results across all batches.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks=2, nlevels=2)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
    "classification", crossval=cval)
  sl <- get_searchlight(ds$dataset, radius=3)
  vox_iter <- lapply(sl, function(x) x)
  results <- mvpa_iterate(mspec, vox_iter[1:5],
    ids=seq_along(vox_iter[1:5]))
#> INFO [2026-04-25 16:08:47] Processing batch 1/1 (5 ROIs in this batch)
#> INFO [2026-04-25 16:08:48] 
#> MVPA Iteration Complete
#> - Total ROIs: 5
#> - Processed: 5
#> - Skipped: 0
# }
```
