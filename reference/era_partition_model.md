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
  vectors are matched to item keys.

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

- ...:

  Additional fields stored on the model spec.

## Value

A model spec of class \`era_partition_model\` for \`run_regional()\` or
\`run_searchlight()\`.
