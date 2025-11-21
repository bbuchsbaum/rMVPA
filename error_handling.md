# Improving Searchlight Error Handling (Proposal)

## Goal
Keep the existing tibble contract (`result`, `indices`, `performance`, `id`, `error`, `error_message`, `warning`, `warning_message`) but standardize how rows are produced and enrich them with a `context` list-column. Downstream code keeps working; upstream failures become transparent and consistent.

## Helpers (new file: `R/result_helpers.R`)
```r
as_ok <- function(val, id, ctx = list()) {
  tibble::tibble(
    result = list(val$result %||% NULL),
    indices = list(val$indices %||% NULL),
    performance = list(val$performance %||% NULL),
    id = id,
    error = FALSE,
    error_message = NA_character_,
    warning = FALSE,
    warning_message = NA_character_,
    context = list(ctx)
  )
}

as_err <- function(msg, id, ctx = list()) {
  tibble::tibble(
    result = list(NULL),
    indices = list(NULL),
    performance = list(NULL),
    id = id,
    error = TRUE,
    error_message = msg,
    warning = FALSE,
    warning_message = NA_character_,
    context = list(ctx)
  )
}

bind_or_err <- function(x, f, id, ctx = list()) {
  if (is.null(x) || inherits(x, "try-error")) {
    msg <- if (inherits(x, "try-error")) as.character(x) else "Unknown error"
    return(as_err(msg, id, ctx))
  }
  tryCatch(f(x), error = function(e) as_err(e$message, id, ctx))
}
```

## Insertion points (low-touch)
1) **`run_future.default` error paths**
   - Replace manual tibble construction for:
     - ROI validation failure (`roi` is NULL).
     - `process_roi` tryCatch error branch.
   - Add context: `list(roi_id = rnum, analysis_type, center_global_id, size = size)`.

2) **`mvpa_iterate` fatal catch**
   - When the outer tryCatch fails, return `as_err("mvpa_iterate fatal error ...", id = NA, ctx = list())` instead of an empty tibble.

3) **Malformed `mvpa_fun` return (randomized/resampled)**
   - Use `as_err` with context: `iter = i`, `has_error_col`, `type`, `raw_str` (the captured `str()`).
   - IDs: `seq_along(cind)`.
   - Leave the good path unchanged.

4) **Optional**: attach `bad_results` as an attribute in standard, randomized, and resampled paths (already done for randomized/resampled) and keep the new `context` list-column; downstream ignores it if unused.

## Why this is minimally disruptive
- Output schema stays the same; combiners and S3 dispatch still work.
- Changes are localized to error-producing sites; no rewrite of processors or combiners.
- Users get structured context to debug “invalid mvpa_fun result” without needing verbose logs.

## Standard context keys
To make debugging predictable, use these keys consistently:
- `roi_id`: ROI/searchlight identifier
- `iter`: Iteration number (for resampling/randomization)
- `analysis_type`: Type of analysis being run
- `center_global_id`: Spatial location of searchlight center
- `size`: ROI size (voxel count)
- `has_error_col`: Whether input had error column
- `type`: Result type (e.g., "randomized", "resampled")
- `raw_str`: Captured structure for malformed results

Sites can add domain-specific keys as needed.

## Rollout steps
1) Add `R/result_helpers.R` with the helpers and `%||%` if needed.
2) Write tests for helper functions (verify schema, context handling).
3) Patch the three sites above to use `as_ok`/`as_err`.
4) Verify existing code expecting exactly 8 columns doesn't break with new `context` column.
5) Keep existing logging; errors now carry structured context even if logs are missed.
6) Update NEWS that failure rows include a `context` list-column for debugging.
