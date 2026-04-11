# Feature-RSA Across States: Domain Adaptation from Encoding to Recall

Feature-RSA is usually demonstrated within one state: you fit a mapping
from a feature space to one ROI, then evaluate that mapping on held-out
data from the same state. Many memory experiments are harder than that.
The feature space is the same across phases, but the brain state
changes. People watch a movie, then recall it. They see an image, then
imagine it. They listen to a story, then retell it.

[`feature_rsa_da_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_da_model.md)
exists for that setting. It lets you keep the feature space fixed,
borrow strength from the source state, and evaluate representational
geometry on held-out target-state folds.

## What does this buy you?

The synthetic example below has the structure you usually want in an
encoding-to-recall analysis:

- domain 1 has clean source-state rows (`X_train`)
- domain 2 has the same feature blocks but a shifted brain mapping
- the target metric is still feature-RSA style geometry on held-out
  target rows

``` r
da_example <- build_da_example()
da_compare <- fit_da_models(da_example)

knitr::kable(da_compare, digits = 3)
```

|             | model       | target_pattern_correlation | target_rdm_correlation | target_r2_full |
|:------------|:------------|---------------------------:|-----------------------:|---------------:|
| coupled_da  | coupled_da  |                      0.963 |                  0.939 |          0.813 |
| builder_da  | builder_da  |                      0.946 |                  0.888 |          0.788 |
| stacked_da  | stacked_da  |                      0.871 |                  0.849 |          0.454 |
| source_only | source_only |                      0.843 |                  0.830 |          0.365 |

In this toy problem, the source-only model already knows something
useful about the feature space, but it misses the state shift. The
adapted models recover more of the target geometry, and the coupled
model does best because it allows the source and target mappings to
differ while still sharing information.

## What lives in each domain?

In a real encoding-to-recall analysis you usually have:

- source rows: the encoding state, often with higher signal or better
  timing
- target rows: the recall state, often noisier and sometimes shorter
- shared feature blocks: the same low-level, semantic, or model-derived
  features in both domains

The design object stores that separation explicitly.

``` r
data.frame(
  source_rows = nrow(da_example$design_fixed$X_train$X),
  target_rows = nrow(da_example$design_fixed$X_test$X),
  feature_columns = ncol(da_example$design_fixed$X_train$X),
  feature_sets = paste(names(da_example$design_fixed$X_train$indices), collapse = ", "),
  target_runs = length(unique(da_example$block_var_test)),
  stringsAsFactors = FALSE
)
#>   source_rows target_rows feature_columns  feature_sets target_runs
#> 1          36          24               6 low, semantic           3
```

[`feature_sets_design()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md)
keeps the source and target predictors separate, and
[`feature_rsa_da_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_da_model.md)
uses held-out target folds for evaluation. If your recall data already
have an externally defined target representation, you can pass it as a
fixed `X_test`.

## How do you fit the adapted model?

The most direct workflow uses a fixed target representation and lets the
model adapt the mapping on target-train rows only.

``` r
coupled_spec <- feature_rsa_da_model(
  dataset = da_example$dataset,
  design = da_example$design_fixed,
  mode = "coupled",
  lambdas = c(low = 0.1, semantic = 0.1),
  alpha_target = 0.3,
  rho = 3,
  rsa_simfun = "spearman"
)

coupled_res <- run_regional(coupled_spec, da_example$region_mask)

coupled_res$performance_table[, c(
  "target_pattern_correlation",
  "target_rdm_correlation",
  "target_r2_full"
)]
#> # A tibble: 1 × 3
#>   target_pattern_correlation target_rdm_correlation target_r2_full
#>                        <dbl>                  <dbl>          <dbl>
#> 1                      0.963                  0.939          0.813
```

Use `mode = "stacked"` when you want one shared mapping estimated from
source rows plus target-train rows. Use `mode = "coupled"` when you
expect a real state shift and want separate source and target weights
tied together by `rho`.

## How do you keep the matching step unbiased?

The important new piece is `target_builder`. Use it when the target-side
representation is not fixed ahead of time, but is itself estimated from
the recall data.

The callback is fold-aware. It receives `train_idx` and `test_idx` from
the outer target fold, and it must return the target predictors in the
original target row order.

``` r
target_builder <- function(X_train, train_idx, builder_data) {
  X_target <- builder_data$X_target
  train_center <- colMeans(X_target[train_idx, , drop = FALSE])
  source_center <- colMeans(X_train$X)
  centered <- sweep(X_target, 2, train_center, "-")
  sweep(0.85 * centered, 2, source_center, "+")
}
```

In this vignette the builder does a simple train-fold-only recentering
step. In a real watch-to-recall pipeline, this is where you would plug
in your row-matching or alignment routine.

``` r
builder_design <- feature_sets_design(
  X_train = da_example$design_fixed$X_train,
  X_test = NULL,
  block_var_test = da_example$block_var_test,
  target_builder = target_builder,
  target_builder_data = list(X_target = da_example$X_rec_noisy),
  n_test = nrow(da_example$X_rec_noisy)
)

builder_spec <- feature_rsa_da_model(
  dataset = da_example$dataset,
  design = builder_design,
  mode = "coupled",
  lambdas = c(low = 0.1, semantic = 0.1),
  alpha_target = 0.3,
  rho = 3,
  rsa_simfun = "spearman"
)

builder_res <- run_regional(builder_spec, da_example$region_mask)

builder_res$performance_table[, c("target_rdm_correlation", "target_r2_full")]
#> # A tibble: 1 × 2
#>   target_rdm_correlation target_r2_full
#>                    <dbl>          <dbl>
#> 1                  0.888          0.788
```

That is the right pattern whenever the target-side alignment would
otherwise leak held-out recall rows into the analysis. The model still
evaluates on held-out target rows, but now the target feature
construction is also fold-aware.

## What changes when recall is a single run?

If the target domain is one continuous run, you still do not have to
abandon the DA workflow. The model will fall back to contiguous target
folds, and you can add a purge gap to reduce temporal leakage.

``` r
single_run_spec <- feature_rsa_da_model(
  dataset = da_example$dataset,
  design = builder_design,
  mode = "coupled",
  lambdas = c(low = 0.1, semantic = 0.1),
  recall_nfolds = 4,
  target_gap = 2,
  target_nperm = 50
)
```

Use that pattern when recall is a single stream rather than a multi-run
design. The contiguous folds define target-train and target-test
segments, `target_gap` removes the rows adjacent to the held-out
segment, and `target_nperm` gives you a single-run null that preserves
local temporal structure.

## When should you use this workflow?

- Use fixed `X_test` when the target representation is external and
  already aligned.
- Use `target_builder` when the target representation is estimated from
  recall itself and must be rebuilt per fold.
- Use `mode = "coupled"` when the source and target mappings should be
  similar but not identical.
- Use `mode = "stacked"` when you want one shared mapping and only mild
  state shift.

Typical use cases include encoding-to-recall transfer,
perception-to-imagery transfer, story listening to retelling, and any
setting where the same feature space survives a change in cognitive
state.

## Next steps

The next question is often no longer about one ROI. If you want to ask
whether feature-RSA geometry generalizes across ROIs, see
[`vignette("Feature_RSA_Connectivity")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Connectivity.md).
If you want the broader decision map for within-ROI fits, cross-state
transfer, and cross-ROI transfer, see
[`vignette("Feature_RSA_Advanced_Workflows")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Advanced_Workflows.md).
