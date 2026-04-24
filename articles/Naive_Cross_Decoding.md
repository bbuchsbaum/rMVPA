# Naive Cross-Decoding

## Overview

You use cross-decoding when the data come from two related domains, such
as perception and memory, and you want to know whether a pattern learned
in the source domain transfers to the target domain. Naive
cross-decoding is the direct-transfer baseline: it classifies
target-domain patterns by correlating them with source-domain
prototypes, with no learned correction between the two domains.

This makes the analysis deliberately conservative. If it works, the
representations are already aligned well enough to transfer directly. If
it does not, the result is still useful because it tells you how much
work a domain-adaptive model such as
[`remap_rrr_model()`](http://bbuchsbaum.github.io/rMVPA/reference/remap_rrr_model.md)
or an encoding-retrieval partition analysis has to explain.

**Algorithm:**

1.  Compute a prototype (mean pattern) for each category in the training
    domain.
2.  Correlate every test pattern with each prototype.
3.  Assign the category whose prototype has the highest correlation.

Correlations are passed through softmax to produce pseudo-probabilities.

## Quick example

``` r
set.seed(42)
data_info <- gen_sample_dataset(
  D = c(10, 10, 10), nobs = 80, nlevels = 4,
  blocks = 4, external_test = TRUE
)

# Three ROIs
roi_mask <- NeuroVol(
  sample(1:3, length(data_info$dataset$mask), replace = TRUE),
  space(data_info$dataset$mask)
)

model <- naive_xdec_model(
  dataset = data_info$dataset,
  design  = data_info$design,
  return_predictions = TRUE
)

print(model)
#> Naive Cross-Decoding Model
#>   link_by:             (class labels) 
#>   training levels:     a, b, c, d 
#>   n_train:             80 
#>   n_test:              80 
#>   return_predictions:  TRUE
```

``` r
results <- run_regional(model, roi_mask)
metric_cols <- intersect(
  c("roinum", "Accuracy", "AUC"),
  names(results$performance_table)
)
results$performance_table[, metric_cols, drop = FALSE]
#> # A tibble: 3 × 3
#>   roinum Accuracy     AUC
#>    <int>    <dbl>   <dbl>
#> 1      1    0.2   -0.0292
#> 2      2    0.238 -0.0562
#> 3      3    0.212  0.0558
```

``` r
chance <- 1 / nlevels(data_info$design$y_test)
cat(sprintf("Chance level: %.0f%%\n", chance * 100))
#> Chance level: 25%
```

With random data, accuracy hovers around chance. With real fMRI data,
above-chance accuracy indicates that category representations transfer
across domains.

## The `link_by` parameter

By default (`link_by = NULL`), prototypes are keyed by class label: all
training trials for category A are averaged into one prototype, and each
test trial is classified against these category prototypes.

When the same individual items appear in both domains (e.g., the same
images during perception and recall), you can form **item-level
prototypes** instead:

``` r
model <- naive_xdec_model(
  dataset = dataset,
  design  = design,
  link_by = "item_id"
)
```

This creates one prototype per item rather than per category, giving
finer-grained matching when item identity is meaningful.

## How does this relate to ERA partitioning?

Naive cross-decoding answers a first-order question: does each target
pattern look most like the correct source prototype? That is often the
first result to check because it maps directly onto classification
accuracy.

The
[`era_partition_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
analysis uses the same direct-transfer idea as one component of a
broader encoding-retrieval analysis. It builds matched source and target
prototypes for the item key, then reports:

- direct transfer metrics such as `naive_top1_acc` and `xdec_Accuracy`,
- first-order variance partitioning via `first_order_delta_r2`,
- second-order geometry preservation via `second_order_delta_r2` and
  `geom_cor`, and
- optional leakage-free Procrustes alignment metrics.

Use
[`naive_xdec_model()`](http://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md)
when you need a clean baseline map or ROI table for direct transfer. Use
[`era_partition_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
when your scientific question is whether encoding-retrieval transfer
reflects item-specific similarity, preserved representational geometry,
or nuisance structure such as block, category, or temporal lag.

## Searchlight analysis

Naive cross-decoding plugs directly into
[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md):

``` r
model <- naive_xdec_model(dataset = dataset, design = design)
sl_results <- run_searchlight(model, radius = 3)
```

The resulting accuracy map shows where category structure transfers
across domains, voxel by voxel.

## Troubleshooting

**“requires external test set”** — Naive cross-decoding needs separate
train and test data. Provide `external_test = TRUE` in
[`gen_sample_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md),
or construct your
[`mvpa_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
with both `train_data` and `test_data`.

**All predictions are the same class** — Prototypes may be too similar.
Check inter-prototype correlations; if all are above ~0.9, the training
domain does not carry enough discriminative information for the method
to work.

**Chance-level performance everywhere** — This is the expected outcome
when there is a large systematic domain shift. Consider domain-adaptive
models such as
[`remap_rrr_model()`](http://bbuchsbaum.github.io/rMVPA/reference/remap_rrr_model.md)
(see
[`vignette("REMAP_RRR")`](http://bbuchsbaum.github.io/rMVPA/articles/REMAP_RRR.md)).

## Next steps

- [`vignette("REMAP_RRR")`](http://bbuchsbaum.github.io/rMVPA/articles/REMAP_RRR.md)
  — domain-adaptive cross-decoding with low-rank corrections
- [`vignette("ERA_RSA_Cross_Decoding")`](http://bbuchsbaum.github.io/rMVPA/articles/ERA_RSA_Cross_Decoding.md)
  — encoding-retrieval RSA and ERA partitioning
- [`?naive_xdec_model`](http://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md)
  — full parameter documentation
