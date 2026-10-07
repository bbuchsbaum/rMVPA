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
# A predictor list as stored in regional_mvpa_result$fits when
# run_regional(..., return_fits = TRUE) is used with remap_rrr_model():
pred <- list(diag_by_fold = list(
  list(train_items = c("a", "b"), item_res_naive = c(4, 8),
       item_res_remap = c(2, 4)),
  list(train_items = c("a", "b"), item_res_naive = c(6, 10),
       item_res_remap = c(3, 5))
))
summarize_remap_items(pred)
#> # A tibble: 2 × 5
#>   item  res_naive res_remap n_folds res_ratio
#>   <chr>     <dbl>     <dbl>   <int>     <dbl>
#> 1 a             5       2.5       2       0.5
#> 2 b             9       4.5       2       0.5
```
