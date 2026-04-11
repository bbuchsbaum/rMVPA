# Naive Cross-Decoding

## Overview

Naive cross-decoding classifies test-domain patterns by correlating them
with training-domain prototypes (class-mean patterns) — no domain
adaptation, no learned transformation. It is the natural baseline for
any cross-domain analysis: if it works, your representations generalise
across contexts without correction.

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
results$performance_table
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
  — encoding–retrieval RSA
- [`?naive_xdec_model`](http://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md)
  — full parameter documentation
