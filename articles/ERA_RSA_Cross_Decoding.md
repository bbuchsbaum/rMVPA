# ERA-RSA: Cross-Decoding Between Encoding and Retrieval

## Overview

The **ERA-RSA model**
([`era_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md))
combines:

- **Encoding–Retrieval Accuracy (ERA)** — first-order similarity between
  encoding and retrieval patterns for each item, and
- **Representational Geometry** — second-order similarity between
  encoding and retrieval representational dissimilarity matrices (RDMs).

This vignette walks through a small synthetic example of **encoding →
retrieval cross-decoding** using ERA-RSA, and illustrates how to run
both **regional** and **searchlight** analyses.

ERA-RSA is designed for situations where you have:

- an encoding phase (e.g., stimulus viewing),
- a retrieval phase (e.g., recognition or recall),
- an item key that links encoding and retrieval trials.

## 1. Synthetic Encoding/Retrieval Dataset

We start from
[`gen_sample_dataset()`](http://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md)
and treat its training split as **encoding** and its external test split
as **retrieval**.

``` r
set.seed(123)

toy <- gen_sample_dataset(
  D            = c(4, 4, 4),  # modest volume for speed
  nobs         = 24,
  nlevels      = 3,           # 3 item categories
  blocks       = 3,           # 3 runs
  external_test = TRUE,       # required for encoding → retrieval
  ntest_obs     = 24
)

# Add an item key column to both train and test designs
item_keys <- sprintf("item_%02d", seq_len(nrow(toy$design$train_design)))
toy$design$train_design$item <- factor(item_keys, levels = item_keys)
toy$design$test_design$item  <- factor(item_keys, levels = item_keys)

head(toy$design$train_design)
#> # A tibble: 6 x 4
#>   Y     block_var .rownum item   
#>   <fct>     <int>   <int> <fct>  
#> 1 b             1       1 item_01
#> 2 a             1       2 item_02
#> 3 a             1       3 item_03
#> 4 a             1       4 item_04
#> 5 c             1       5 item_05
#> 6 b             1       6 item_06
```

Each row in `train_design` and `test_design` corresponds to a trial. The
`item` column links encoding and retrieval presentations of the same
stimulus. The toy data use paired row IDs as item keys so the example
has one encoding and one retrieval observation per item; in an applied
analysis, this column would usually be a stimulus, image, word, or event
identifier from your design table.

### `key_var` vs `y_train` / `y_test` and `phase_var`

An `mvpa_design` already contains `y_train` / `y_test`, which are the
**response labels for each trial** (often category labels). ERA-RSA adds
two extra notions:

- **`key_var`**: an *item key* that links encoding and retrieval trials
  belonging to the same underlying stimulus (e.g., image ID). ERA
  metrics such as `era_top1_acc` operate at the item level: for each
  retrieval trial we ask whether its top match among encoding trials has
  the same `key_var`, not just the same category.
- **`phase_var`**: a phase label (e.g., encoding vs retrieval). In the
  external-test setup used here, the train split is treated as encoding
  and the test split as retrieval, so `phase_var` is mostly a
  placeholder to keep the interface consistent. It becomes important in
  single-dataset designs where both phases live in the same time series.

In short, `y_train`/`y_test` describe *what was shown* on each trial
(class/category), while `key_var` defines *which item* that trial
belongs to; ERA-RSA uses the item keys to build encoding–retrieval
similarity and geometry.

## 2. Regional ERA-RSA Model

We build an `era_rsa_model` and run a simple regional analysis with a
small number of regions.

The two similarity arguments operate at different levels:

- `distfun = cordist("pearson")` says how to build each within-phase
  RDM. Within a region, ERA-RSA first averages trials by `key_var`, then
  computes one encoding RDM and one retrieval RDM from those item-level
  patterns.
- `rsa_simfun = "spearman"` says how to compare those two RDMs. This is
  the second-order RSA step: the lower triangle of the encoding RDM is
  correlated with the lower triangle of the retrieval RDM. The resulting
  regional metric is `geom_cor`.

There is no external model RDM in this example. The “model” for the
second-order comparison is the encoding geometry itself: we ask whether
the relative distances among encoded items are preserved during
retrieval.

``` r
# Simple 3-region mask for demonstration
region_mask <- NeuroVol(
  sample(1:3, size = length(toy$dataset$mask), replace = TRUE),
  space(toy$dataset$mask)
)

era_ms <- era_rsa_model(
  dataset   = toy$dataset,
  design    = toy$design,
  key_var   = ~ item,          # item key linking encoding ↔ retrieval
  phase_var = ~ block_var,     # parsed for interface; train/test define phases here
  distfun   = cordist("pearson"), # builds encoding and retrieval RDMs
  rsa_simfun = "spearman"         # compares those RDMs; output is geom_cor
)

era_ms
#> 
#>  Model Specification 
#> 
#> - Summary 
#>   - Name:  era_rsa_model 
#>   - Primary Class:  era_rsa_model 
#>   - Class Chain:  era_rsa_model, model_spec, list 
#>   - Return Predictions:  FALSE 
#>   - Compute Performance:  TRUE 
#>   - Has Test Set:  TRUE 
#>   - Basis Count:  1 
#> - Parameters 
#>   - key_var: ~item
#>   - phase_var: ~block_var
#>   - encoding_level: 1
#>   - retrieval_level: 2
#>   - distfun: cordist / pearson
#>   - rsa_simfun: spearman
#>   - include_diag: TRUE
#> 
#>  MVPA Dataset 
#> 
#> - Training Data 
#>   - Dimensions:  4 x 4 x 4 x 24 observations 
#>   - Type:  DenseNeuroVec 
#> - Test Data 
#>   - Dimensions:  4 x 4 x 4 x 24 observations 
#>   - Type:  DenseNeuroVec 
#> - Mask Information 
#>   - Areas:  TRUE : 64 
#>   - Active voxels/vertices:  64 
#> 
#> 
#>  MVPA Design 
#> 
#> - Training Data 
#>   - Observations:  24 
#>   - Response Type:  Factor
#>   - Levels:  a, b, c 
#>   - Class Distribution:  a: 8, b: 8, c: 8 
#> - Test Data 
#>   - Observations:  24 
#>   - Class Distribution:  a: 8, b: 8, c: 8 
#> - Structure 
#>   - Blocking:  Present
#>   - Number of Blocks:  3 
#>   - Mean Block Size:  8  (SD:  0 ) 
#>   - Split Groups:  None

era_res <- run_regional(era_ms, region_mask)
era_res
#> 
#>  Regional Analysis Results 
#> 
#> - Summary 
#>   - Model:  era_rsa_model 
#>   - Regions with Results:  3 
#>   - Metrics:  n_items, era_top1_acc, era_diag_mean, era_diag_minus_off, geom_cor, era_diag_minus_off_same_block ... 
#>   - Output Maps:  n_items, era_top1_acc, era_diag_mean, era_diag_minus_off, geom_cor, era_diag_minus_off_same_block ...
```

The `performance_table` contains one row per region and several ERA-RSA
metrics.

``` r
head(era_res$performance_table)
#> # A tibble: 3 x 11
#>   roinum n_items era_top1_acc era_diag_mean era_diag_minus_off geom_cor
#>    <int>   <dbl>        <dbl>         <dbl>              <dbl>    <dbl>
#> 1      1      24       0.0417       -0.0153            -0.0172  -0.0539
#> 2      2      24       0.0833       -0.0437            -0.0481   0.0179
#> 3      3      24       0.0417        0.0509             0.0387   0.0304
#> # i 5 more variables: era_diag_minus_off_same_block <dbl>,
#> #   era_diag_minus_off_diff_block <dbl>, era_lag_cor <dbl>,
#> #   geom_cor_run_partial <dbl>, geom_cor_xrun <dbl>
```

Key metrics include:

- `era_top1_acc` — top-1 encoding→retrieval accuracy for the item key.
- `era_diag_mean` — mean encoding–retrieval similarity on the diagonal.
- `era_diag_minus_off` — diagonal minus off-diagonal similarity.
- `geom_cor` — the second-order correlation between the encoding and
  retrieval RDMs. With `rsa_simfun = "spearman"`, this is a rank
  correlation over item-pair dissimilarities.

These quantify both cross-decoding performance and representational
geometry alignment between phases.

``` r
era_res$performance_table[
  ,
  c("roinum", "era_top1_acc", "era_diag_minus_off", "geom_cor")
]
#> # A tibble: 3 x 4
#>   roinum era_top1_acc era_diag_minus_off geom_cor
#>    <int>        <dbl>              <dbl>    <dbl>
#> 1      1       0.0417            -0.0172  -0.0539
#> 2      2       0.0833            -0.0481   0.0179
#> 3      3       0.0417             0.0387   0.0304
```

## 3. When do you need ERA partitioning?

The ERA-RSA summary above gives you direct item matching and a raw
geometry correlation.
[`era_partition_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
asks a sharper follow-up question: how much of the encoding-retrieval
relationship is uniquely explained by same-item transfer, and how much
is uniquely explained by preserved representational geometry after
nuisance structure is accounted for?

This makes it useful when a single cross-decoding accuracy is too
compressed. For example, two regions could have similar top-1 accuracy,
but one may carry item-specific diagonal similarity while another mainly
preserves the relative geometry among items.

``` r
partition_ms <- era_partition_model(
  dataset = toy$dataset,
  design  = toy$design,
  key_var = ~ item,
  distfun = cordist("pearson"),
  rsa_simfun = "spearman",
  xdec_link_by = "item",
  include_procrustes = FALSE
)

partition_res <- run_regional(partition_ms, region_mask)
```

``` r
partition_cols <- c(
  "roinum",
  "n_items",
  "naive_top1_acc",
  "xdec_Accuracy",
  "first_order_delta_r2",
  "second_order_delta_r2",
  "geom_cor"
)
partition_res$performance_table[, partition_cols]
#> # A tibble: 3 x 7
#>   roinum n_items naive_top1_acc xdec_Accuracy first_order_delta_r2
#>    <int>   <dbl>          <dbl>         <dbl>                <dbl>
#> 1      1      24         0.0417        0.0417             0.000244
#> 2      2      24         0.0833        0.0833             0.00198 
#> 3      3      24         0.0417        0.0417             0.00123 
#> # i 2 more variables: second_order_delta_r2 <dbl>, geom_cor <dbl>
```

Read these metrics in layers:

- `naive_top1_acc` is the item-level direct-transfer score computed from
  the encoding-by-retrieval similarity matrix.
- `xdec_Accuracy` is the trial-level direct cross-decoding metric,
  matched to the same prototype scorer used by
  [`naive_xdec_model()`](http://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md).
- `first_order_delta_r2` is the unique variance explained by same-item
  transfer in the cross-state similarity matrix.
- `second_order_delta_r2` and `geom_cor` summarize whether the encoding
  item geometry is preserved at retrieval.

For this random toy dataset, the values should be near chance or near
zero. In real data, a region with high `first_order_delta_r2` and weak
`second_order_delta_r2` supports item-specific reinstatement without
much global geometry preservation. The reverse pattern suggests that the
relative arrangement of items is preserved even when direct item
identification is weak.

[`era_partition_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
can also include block, category, temporal, and custom nuisance
regressors through `item_block_enc`, `item_block_ret`, `item_time_enc`,
`item_time_ret`, `item_category`, `first_order_nuisance`, and
`second_order_nuisance`. If you have enough paired items, set
`include_procrustes = TRUE` to add a leakage-free orthogonal alignment
decoder.

`include_procrustes` controls an optional decoder that first learns an
orthogonal Procrustes map from encoding prototypes into retrieval
prototype space. The fit is leave-one-item-out: when scoring item *i*,
the rotation is estimated from all other paired items, then the held-out
retrieval item is matched against the aligned encoding prototypes. This
gives Procrustes metrics such as `procrustes_top1_acc`,
`procrustes_diag_mean`, and `procrustes_diag_minus_off`. Use it when you
want to ask whether encoding and retrieval occupy the same
representational space up to a global rotation or reflection. Leave it
off for small toy examples or when there are too few paired items for
stable held-out alignment; `min_procrustes_train_items` sets that
minimum.

## 4. Searchlight ERA-RSA

We can also run ERA-RSA in a searchlight mode to obtain whole-brain maps
of the same metrics.

``` r
set.seed(456)

era_sl <- run_searchlight(
  era_ms,
  radius = 2,         # searchlight radius in voxels
  method = "standard" # ERA-RSA currently uses standard searchlight
)

era_sl
#> 
#>  Searchlight Analysis Results 
#> 
#> - Coverage 
#>   - Voxels/Vertices in Mask:  64 
#>   - Voxels/Vertices with Results:  64 
#> - Output Maps (Metrics) 
#>   -  n_items  (Type:  DenseNeuroVol ) 
#>   -  era_top1_acc  (Type:  DenseNeuroVol ) 
#>   -  era_diag_mean  (Type:  DenseNeuroVol ) 
#>   -  era_diag_minus_off  (Type:  DenseNeuroVol ) 
#>   -  geom_cor  (Type:  DenseNeuroVol ) 
#>   -  era_diag_minus_off_same_block  (Type:  DenseNeuroVol ) 
#>   -  era_diag_minus_off_diff_block  (Type:  DenseNeuroVol ) 
#>   -  era_lag_cor  (Type:  DenseNeuroVol ) 
#>   -  geom_cor_run_partial  (Type:  DenseNeuroVol ) 
#>   -  geom_cor_xrun  (Type:  DenseNeuroVol )
```

The `searchlight_result` contains:

- `metrics`: names of the output maps (e.g., `geom_cor`,
  `era_top1_acc`),
- `results`: a list of `NeuroVol` maps, one per metric.

``` r
era_sl$metrics
#>  [1] "n_items"                       "era_top1_acc"                 
#>  [3] "era_diag_mean"                 "era_diag_minus_off"           
#>  [5] "geom_cor"                      "era_diag_minus_off_same_block"
#>  [7] "era_diag_minus_off_diff_block" "era_lag_cor"                  
#>  [9] "geom_cor_run_partial"          "geom_cor_xrun"
```

We can save the searchlight maps using
[`save_results()`](http://bbuchsbaum.github.io/rMVPA/reference/save_results.md):

``` r
out_dir <- tempfile("era_rsa_sl_")
dir.create(out_dir, showWarnings = FALSE)

save_results(era_sl, out_dir, level = "standard")
list.files(file.path(out_dir, "maps"))
```

This will create one NIfTI file per metric (e.g., `geom_cor.nii.gz`,
`era_top1_acc.nii.gz`) that can be viewed in your favorite neuroimaging
software.

## 5. Adding Confounds and Lag Information

ERA-RSA can optionally incorporate item-level confounds and lag
variables. These are all defined at the **item** level, not the trial
level, and must align with the levels of `key_var`.

### Item-level confounds (`confound_rdms`)

`confound_rdms` is a named list of K×K matrices or `"dist"` objects
describing item-by-item nuisance structure (e.g., block/run/time), where
rows/columns correspond to item keys.

A common pattern is to build run-based confounds:

``` r
items <- levels(toy$design$train_design$item)

# Example: per-item encoding run (modal run for each item)
Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
item_run_enc <- sapply(items, function(it) {
  Mode(toy$design$train_design$block_var[toy$design$train_design$item == it])
})
names(item_run_enc) <- items

# Example: per-item retrieval run (could differ from encoding)
item_run_ret <- sample(item_run_enc)

# Build run confound RDMs on the item grid
run_enc <- outer(item_run_enc, item_run_enc, FUN = function(a, b) as.numeric(a == b))
run_ret <- outer(item_run_ret, item_run_ret, FUN = function(a, b) as.numeric(a == b))
rownames(run_enc) <- colnames(run_enc) <- items
rownames(run_ret) <- colnames(run_ret) <- items

era_ms_conf <- era_rsa_model(
  dataset = toy$dataset,
  design  = toy$design,
  key_var = ~ item,
  phase_var = ~ block_var,
  confound_rdms = list(run_enc = run_enc, run_ret = run_ret),
  item_run_enc  = factor(item_run_enc),
  item_run_ret  = factor(item_run_ret)
)
```

With these supplied, two additional metrics become available:

- `geom_cor_run_partial`: correlation between encoding and retrieval
  RDMs after regressing out the enc/ret run RDMs.
- `geom_cor_xrun`: correlation between encoding and retrieval RDMs
  restricted to item pairs that differ in both encoding and retrieval
  run.

### Block structure (`item_block`)

`item_block` encodes a per-item block/condition, typically derived from
`design$train_design$block_var`:

``` r
item_block <- factor(item_run_enc, levels = sort(unique(item_run_enc)))

era_ms_block <- era_rsa_model(
  dataset    = toy$dataset,
  design     = toy$design,
  key_var    = ~ item,
  phase_var  = ~ block_var,
  item_block = item_block
)
```

This enables block-limited ERA contrasts:

- `era_diag_minus_off_same_block`: diagonal ERA minus mean similarity to
  other items in the same block.
- `era_diag_minus_off_diff_block`: diagonal ERA minus mean similarity to
  items in different blocks.

### Lag information (`item_lag`)

`item_lag` is a numeric per-item lag between encoding and retrieval
(e.g., retrieval onset – encoding onset), aligned to item keys. It is
used to compute:

- `era_lag_cor`: Spearman correlation between the ERA diagonal and
  `item_lag`, using complete cases.

In applied analyses, you would construct `item_block`, `item_lag`,
`item_run_enc`, `item_run_ret`, and `confound_rdms` from your
experiment’s design tables following the patterns above, then pass them
into
[`era_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md).

## 6. Summary

- [`era_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md)
  provides a unified framework for:
  - cross-decoding between encoding and retrieval, and
  - comparing encoding and retrieval representational geometries.
- [`era_partition_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
  separates direct item transfer, variance uniquely explained by
  same-item similarity, and variance explained by preserved second-order
  geometry.
- It integrates naturally with:
  - [`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
    for ROI-based analyses, and
  - [`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
    for whole-brain mapping.
- Outputs are standard rMVPA result objects, so the same tooling
  ([`save_results()`](http://bbuchsbaum.github.io/rMVPA/reference/save_results.md),
  plotting, etc.) applies as for other models such as
  [`naive_xdec_model()`](http://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md)
  and
  [`remap_rrr_model()`](http://bbuchsbaum.github.io/rMVPA/reference/remap_rrr_model.md).

ERA-RSA is particularly useful when you want to go beyond simple
accuracy and ask **how similar the geometry of neural representations is
across phases**, and how that similarity varies across brain regions.
