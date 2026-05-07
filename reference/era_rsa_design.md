# Build item-level ERA-RSA confounds and summaries from mvpa_design

Constructs per-item summaries (block/run/time/lag) and item-level
confound RDMs that plug directly into
`era_rsa_model(..., confound_rdms=...)`. Reuses temporal helpers in
R/temporal_rdms.R to avoid redundancy.

## Usage

``` r
era_rsa_design(
  design,
  key_var,
  phase_var,
  encoding_level = NULL,
  retrieval_level = NULL,
  block_var = NULL,
  time_var = NULL,
  phase_scoped_runs = FALSE
)
```

## Arguments

- design:

  mvpa_design (same object passed to era_rsa_model).

- key_var:

  Column name or formula giving the item ID (e.g., ~ ImageID).

- phase_var:

  Column name or formula indicating phase (e.g., ~ Phase).

- encoding_level:

  Encoding phase label; default = first level of phase_var.

- retrieval_level:

  Retrieval phase label; default = second level of phase_var.

- block_var:

  Optional column giving block/run membership.

- time_var:

  Optional column giving trial index or onset time.

- phase_scoped_runs:

  Logical; when `TRUE`, prefix per-item run labels with `enc_` / `ret_`
  so encoding and retrieval scans with overlapping numeric labels (e.g.
  both phases have runs 1, 2, 3) are not treated as the same run during
  cross-phase ERA comparisons. Default `FALSE` for back-compatibility;
  set to `TRUE` when the encoding and retrieval phases come from
  different scan sessions.

## Value

A list with elements:

- items: factor of item IDs used (intersection of E/R)

- item_block: optional factor of per-item block labels

- item_time_enc, item_time_ret: optional numeric vectors

- item_lag: optional numeric vector (ret - enc)

- item_run_enc, item_run_ret: optional per-item run factors

- confound_rdms: named list of item-level RDMs (block/time/run)

## Examples

``` r
if (FALSE) { # \dontrun{
  # See era_rsa_model for full ERA-RSA workflow with design construction
} # }
```
