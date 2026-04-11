# Summarize REMAP item-level residuals for an ROI

Extract per-item distances between memory and perception prototypes in
the jointly-whitened space, before (naive) and after the REMAP
correction. Works from a \`regional_mvpa_result\` when
\`return_fits=TRUE\`, or directly from a predictor list returned by
\`process_roi()\`.

## Usage

``` r
summarize_remap_items(x, roi = NULL)
```

## Arguments

- x:

  Either a \`regional_mvpa_result\` (preferred) or a predictor list
  containing \`diag_by_fold\`.

- roi:

  ROI index or id when \`x\` is a regional result. If \`NULL\` and \`x\`
  is a predictor, it is ignored.

## Value

A tibble with columns: \`item\`, \`res_naive\`, \`res_remap\`,
\`res_ratio\`, and \`n_folds\`.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- run_regional(model, region_mask, return_fits = TRUE)
items_tbl <- summarize_remap_items(res, roi = 1)
} # }
```
