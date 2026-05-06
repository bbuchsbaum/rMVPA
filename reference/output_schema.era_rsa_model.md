# Output schema for era_rsa_model

Declares all scalar metrics emitted by `fit_roi.era_rsa_model` so that
`combine_schema_standard` can pre-allocate output maps with consistent
length across ROIs.

## Usage

``` r
# S3 method for class 'era_rsa_model'
output_schema(model)
```

## Arguments

- model:

  An era_rsa_model object.

## Value

A named character vector mapping each emitted metric name to `"scalar"`.
Base metrics are `n_items`, `era_top1_acc`, `era_diag_mean`,
`era_diag_minus_off`, `geom_cor`, `era_diag_minus_off_same_block`,
`era_diag_minus_off_diff_block`, `era_lag_cor`, `geom_cor_run_partial`,
and `geom_cor_xrun`. If `confound_rdms` is supplied, the schema also
includes `beta_enc_geom`, one `beta_<name>` per confound RDM,
`sp_enc_geom`, and one `sp_<name>` per confound RDM.
