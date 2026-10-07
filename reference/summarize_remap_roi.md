# Summarize REMAP diagnostics at the ROI level

Convenience helper to extract ROI-level REMAP diagnostics from a
\`regional_mvpa_result\`. It prefers metrics recorded in the
\`performance_table\` (so it works even when fits are not returned), and
will fall back to averaging diagnostics stored in
\`fits\[\[i\]\]\$diag_by_fold\` when needed.

## Usage

``` r
summarize_remap_roi(regional_res)
```

## Arguments

- regional_res:

  A \`regional_mvpa_result\` returned by \`run_regional()\`.

## Value

A tibble with one row per ROI containing: - \`roinum\`: ROI id -
\`mean_rank\`: mean adapter rank - \`mean_lambda\`: mean selected
lambda - \`mean_roi_improv\`: mean fraction of P-\>M mismatch
explained - \`mean_delta_frob\`: mean Frobenius norm of the learned
correction

## Examples

``` r
# Normally `res` comes from run_regional() on a remap_rrr_model(); here a
# minimal regional result carrying REMAP columns in its performance table:
res <- structure(
  list(performance_table = data.frame(
    roinum = 1:2, adapter_rank = c(2, 3), lambda_mean = c(0.25, 0.5),
    remap_improv = c(0.30, 0.45), delta_frob_mean = c(1.2, 1.8)
  )),
  class = c("regional_mvpa_result", "list")
)
summarize_remap_roi(res)
#> # A tibble: 2 × 5
#>   roinum mean_rank mean_lambda mean_roi_improv mean_delta_frob
#>    <int>     <dbl>       <dbl>           <dbl>           <dbl>
#> 1      1         2        0.25            0.3              1.2
#> 2      2         3        0.5             0.45             1.8
```
