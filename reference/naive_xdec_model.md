# Naive Cross-Decoding (correlation to source prototypes)

Baseline analysis that classifies target-domain trials by correlating
them directly with source-domain prototypes (means over repeats) without
any adaptation. Serves as a comparator for REMAP-RRR.

## Usage

``` r
naive_xdec_model(
  dataset,
  design,
  link_by = NULL,
  return_predictions = TRUE,
  ...
)
```

## Arguments

- dataset:

  mvpa_dataset with train_data (source) and test_data (target)

- design:

  mvpa_design with y_train/y_test and optional \`link_by\` column
  present in both train/test designs. If \`link_by\` is NULL, prototypes
  are formed by the training labels \`y_train\` and evaluated against
  \`y_test\`.

- link_by:

  Optional character; if provided, prototypes and observed labels are
  keyed by this column instead of y.

- return_predictions:

  logical; keep per-ROI predictions.

- ...:

  Additional arguments (currently unused).

## Value

model spec object of class \`naive_xdec_model\` for use with
\`run_regional()\` or \`run_searchlight()\`.

## Details

The term "naive" is a technical designation meaning "without domain
adaptation" (i.e., direct prototype matching). It serves as the
principled baseline for evaluating domain-adaptive methods like
REMAP-RRR. Also known as "zero-shot transfer" or "direct transfer" in
machine learning literature.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires dataset with train_data and test_data
  ds <- gen_sample_dataset(c(5,5,5), 20, external_test=TRUE)
  model <- naive_xdec_model(ds$dataset, ds$design)
} # }
```
