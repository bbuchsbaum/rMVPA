# Run Future

Run a future-based computation defined by the object and frame.

## Usage

``` r
run_future(obj, frame, processor, ...)

# Default S3 method
run_future(
  obj,
  frame,
  processor = NULL,
  verbose = FALSE,
  analysis_type = "searchlight",
  drop_probs = FALSE,
  fail_fast = FALSE,
  item_seeds = NULL,
  ...
)

# S3 method for class 'shard_model_spec'
run_future(
  obj,
  frame,
  processor = NULL,
  verbose = FALSE,
  analysis_type = "searchlight",
  drop_probs = FALSE,
  fail_fast = FALSE,
  item_seeds = NULL,
  ...
)
```

## Arguments

- obj:

  An object specifying the computation.

- frame:

  A data frame or environment containing data for the computation.

- processor:

  A function or object specifying how to process the frame.

- ...:

  Additional arguments passed to the method-specific function.

- verbose:

  Logical; print progress messages if `TRUE`.

- analysis_type:

  The type of analysis (e.g., "searchlight").

- drop_probs:

  Logical; if `TRUE`, drop per-observation class probabilities from
  returned results and keep only the observed-class probability.

- fail_fast:

  Logical; if `TRUE`, stop on the first ROI error.

- item_seeds:

  Optional list with one L'Ecuyer-CMRG RNG stream (`.Random.seed`
  vector) per row of `frame`.
  [`mvpa_iterate()`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
  supplies the streams it drew for the whole iteration, so a ROI's
  random numbers do not depend on batching or the future plan. When
  `NULL`, streams are drawn here.

## Value

A tibble containing the processed results.

## Details

If the progressr package is installed, this method emits per-task
progress updates from parallel workers. Progress handling is enabled for
the scope of this call and uses the currently configured handlers.

## Examples

``` r
frame <- tibble::tibble(
  .id = 1:2,
  rnum = c("roi1", "roi2"),
  roi = list(1:3, 4:5),
  size = c(3, 2)
)
mod_spec <- list(process_roi = function(mod_spec, roi, rnum, ...) {
  tibble::tibble(
    result = list(mean(roi)),
    indices = list(seq_along(roi)),
    performance = list(NULL),
    id = rnum
  )
})
run_future(mod_spec, frame, NULL)
#> Warning: process_roi legacy dispatch is deprecated. Implement a fit_roi.* method for your model class.
#> Warning: process_roi legacy dispatch is deprecated. Implement a fit_roi.* method for your model class.
#> # A tibble: 2 × 8
#>   result indices performance id    error error_message   warning warning_message
#>   <list> <list>  <list>      <chr> <lgl> <chr>           <lgl>   <chr>          
#> 1 <NULL> <NULL>  <NULL>      roi1  TRUE  "Error process… TRUE    Error processi…
#> 2 <NULL> <NULL>  <NULL>      roi2  TRUE  "Error process… TRUE    Error processi…
```
