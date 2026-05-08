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
  performance = NULL,
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

- performance:

  Optional user-supplied performance function. When non-\`NULL\`, must
  be a function that accepts the classification result object (with
  fields \`observed\`, \`predicted\`, \`probs\`, \`testind\`,
  \`test_design\`) and returns a named numeric vector of metrics. Routed
  through \[get_custom_perf()\], which annotates it as \`"custom"\`;
  this disables the optimised fast-metric kernel and ensures the user's
  function is always called on the full result object. When \`NULL\`,
  the default binary / multiclass / regression performance helper is
  used.

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

  # Custom metric that uses a column from test_design
  custom_fun <- function(result) {
    vivid <- result$test_design$RateVivid
    probs <- as.matrix(result$probs)
    obs   <- as.character(result$observed)
    true_p <- probs[cbind(seq_along(obs), match(obs, colnames(probs)))]
    c(vivid_spearman = stats::cor(vivid, true_p, method = "spearman"))
  }
  model <- naive_xdec_model(ds$dataset, ds$design, performance = custom_fun)
} # }
```
