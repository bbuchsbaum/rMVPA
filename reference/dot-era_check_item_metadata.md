# Validate item-level metadata for ERA models

Used by
[`era_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md)
and
[`era_partition_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
to warn (or error, when `strict=TRUE`) about missing item-level vectors
that would otherwise produce silently-NA metrics or unused nuisance
regressors. Also detects an easy-to-miss namespace collision where
encoding and retrieval run labels share atomic values without being
phase-scoped.

## Usage

``` r
.era_check_item_metadata(
  where,
  item_block = NULL,
  item_run_enc = NULL,
  item_run_ret = NULL,
  strict = FALSE,
  metric_for_block = character(),
  metric_for_run = character()
)
```

## Arguments

- where:

  Character; calling function name used in messages.

- item_block, item_run_enc, item_run_ret:

  Optional item-level vectors.

- strict:

  Logical; if `TRUE`, missing metadata becomes an error.

- metric_for_block:

  Character vector of metric names that depend on `item_block`.

- metric_for_run:

  Character vector of metric names that depend on `item_run_enc` /
  `item_run_ret`.

## Value

Invisibly `NULL`; called for side effects.
