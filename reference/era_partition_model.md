# ERA Variance-Partition Model

Compares first-order encoding-retrieval transfer with second-order
geometry preservation using matched variance-partition models. The model
builds item-level prototypes in a source state (\`dataset\$train_data\`)
and a target state (\`dataset\$test_data\`), then estimates how much
variance is uniquely explained by same-item transfer and by source-state
geometry after optional nuisance pair models.

## Usage

``` r
era_partition_model(
  dataset,
  design,
  key_var,
  distfun = cordist(method = "pearson"),
  rsa_simfun = c("pearson", "spearman"),
  first_order_nuisance = NULL,
  second_order_nuisance = NULL,
  item_block_enc = NULL,
  item_block_ret = NULL,
  item_run_enc = NULL,
  item_run_ret = NULL,
  item_time_enc = NULL,
  item_time_ret = NULL,
  item_category = NULL,
  compute_xdec_performance = TRUE,
  xdec_link_by = NULL,
  include_procrustes = TRUE,
  procrustes_center = TRUE,
  min_procrustes_train_items = 3L,
  return_matrices = FALSE,
  return_xdec_predictions = FALSE,
  auto_nuisance = TRUE,
  global_nuisance = FALSE,
  require_run_metadata = FALSE,
  ...
)
```

## Arguments

- dataset:

  An \`mvpa_dataset\` with \`train_data\` and \`test_data\`.

- design:

  An \`mvpa_design\` with train/test design tables.

- key_var:

  Column name or formula giving the item key shared across source and
  target states.

- distfun:

  Distance function used for within-state RDMs.

- rsa_simfun:

  Correlation method for the raw geometry summary.

- first_order_nuisance:

  Optional named list of \`K x K\` matrices or length \`K^2\` vectors
  for cross-state similarity nuisance regressors. Matrices are
  interpreted as target rows by source columns.

- second_order_nuisance:

  Optional named list of \`K x K\` matrices or lower-triangle vectors
  for geometry nuisance regressors.

- item_block_enc, item_block_ret:

  Optional item-level block labels for source and target states. Named
  vectors are matched to item keys. Both are required to enable the
  same/different block nuisance regressors used by
  `.era_partition_first_nuisance()` and
  `.era_partition_second_nuisance()`; supplying only one is treated as
  missing for cross-state nuisance purposes.

- item_run_enc, item_run_ret:

  Optional item-level run labels for source and target states. Named
  vectors are matched to item keys and, when `auto_nuisance` includes
  `"run"`, add `same_run_cross`, `same_run_enc`, and `same_run_ret`
  nuisance regressors.

- item_time_enc, item_time_ret:

  Optional item-level time/order values for source and target states.
  Named vectors are matched to item keys.

- item_category:

  Optional item-level category labels used to add a same-category
  nuisance model to both first- and second-order regressions.

- compute_xdec_performance:

  Logical; compute trial-level naive cross-decoding performance using
  the same prototype scorer as
  [`naive_xdec_model`](http://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md).

- xdec_link_by:

  Optional column name used to define source/target labels for the
  trial-level cross-decoding metrics. If `NULL`, `key_var` is used.

- include_procrustes:

  Logical; compute leave-one-item-out orthogonal Procrustes
  cross-decoding metrics.

- procrustes_center:

  Logical; center source and target prototypes using only
  alignment-training items before fitting Procrustes maps.

- min_procrustes_train_items:

  Minimum number of paired items allowed for each leave-one-item-out
  Procrustes alignment.

- return_matrices:

  Logical; store prototype/similarity matrices in each ROI result for
  diagnostics.

- return_xdec_predictions:

  Logical; store the trial-level `classification_result` produced by the
  direct cross-decoder in each ROI result.

- auto_nuisance:

  Logical or character vector controlling automatically derived
  item-level nuisance regressors. `TRUE` includes available `"block"`,
  `"run"`, `"time"`, `"category"`, and `"global"` regressors. `FALSE`
  disables all automatic nuisance regressors so only
  `first_order_nuisance` and `second_order_nuisance` are used. A
  character vector selects specific groups.

- global_nuisance:

  Logical or pre-supplied list controlling whole-mask global similarity
  nuisance. `FALSE` (default) disables it. `TRUE` computes item-level
  whole-mask similarity/RDMs over `dataset$mask` once at construction
  time. A pre-computed list with `S_cross`/ `first`, `D_enc`/`enc`, and
  `D_ret`/`ret` matrices can be supplied directly. When `auto_nuisance`
  includes `"global"`, the cross-state similarity enters the first-order
  model as `global_cross`, and the encoding/retrieval RDMs enter the
  second-order model as `global_enc` and `global_ret`. Caveat: each
  ROI/sphere is part of the global mask, so for large regional ROIs
  covering most of the mask the residualization partially removes local
  signal too.

- require_run_metadata:

  Logical; if `TRUE`, missing item-level block metadata becomes an error
  rather than a warning. Use this when downstream nuisance partitioning
  depends on the same/different-block regressors. Default `FALSE` (warn
  only).

- ...:

  Additional fields stored on the model spec.

## Value

A model spec of class \`era_partition_model\` for \`run_regional()\` or
\`run_searchlight()\`.

## Trial-level vs. item-level metadata

`block_var` on
[`mvpa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
is *trial-level* and is not automatically used here. The first- and
second-order nuisance regressors (`same_block_cross`, `same_block_enc`,
`same_block_ret`) require *item-level* vectors named by levels of
`key_var`: pass `item_block_enc` and `item_block_ret` explicitly. Run
nuisance regressors similarly require `item_run_enc` and `item_run_ret`.
Use `auto_nuisance = FALSE` to suppress all automatic
block/run/time/category regressors when supplying a custom nuisance
model, or pass a character vector such as `c("run", "time")` to keep
only selected groups. When run labels overlap across encoding and
retrieval scans, use phase-scoped labels such as `enc_1` / `ret_1` so
item-level equality across phases is meaningful.
