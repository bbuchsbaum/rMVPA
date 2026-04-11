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
if (FALSE) { # \dontrun{
res <- run_regional(model, region_mask, return_fits = TRUE)
summarize_remap_roi(res)
} # }
```
