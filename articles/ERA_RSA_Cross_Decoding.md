# ERA-RSA: Cross-Decoding Between Encoding and Retrieval

## Overview

The **ERA-RSA model**
([`era_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md))
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
[`gen_sample_dataset()`](https://bbuchsbaum.github.io/rMVPA/reference/gen_sample_dataset.md)
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
  and the test split as retrieval, so `phase_var` is optional. It is
  required when both phases live in the same image series.

In short, `y_train`/`y_test` describe *what was shown* on each trial
(class/category), while `key_var` defines *which item* that trial
belongs to; ERA-RSA uses the item keys to build encoding–retrieval
similarity and geometry.

For a combined series with alternating encoding and retrieval blocks,
keep one design row per image volume and identify phases explicitly:

``` r

era_ms <- era_rsa_model(
  dataset = combined_dataset,
  design = combined_design,
  key_var = ~ item,
  phase_var = ~ phase,
  encoding_level = "E",
  retrieval_level = "R",
  pairing = "one_to_one"
)
```

Pairing is always by `key_var`, never by row adjacency or block order.
rMVPA splits the image series internally while preserving the original
observation order within each phase.

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
  pairing   = "one_to_one",    # require one encoding and one retrieval row per item
  era_simfun = "pearson",      # item-level encoding-retrieval similarity
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
#>   - encoding_level: encoding
#>   - retrieval_level: retrieval
#>   - pairing: one_to_one
#>   - era_simfun: pearson
#>   - era_min_voxels: 2
#>   - era_components: item, identification, geometry
#>   - era_association_score: item_specificity
#>   - era_cor_method: spearman
#>   - era_min_complete: 4
#>   - era_combined_input: FALSE
#>   - distfun: cordist / pearson
#>   - rsa_simfun: spearman
#>   - partial_against: run
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
#> # A tibble: 3 x 12
#>   roinum n_items era_top1_acc era_diag_mean era_diag_minus_off geom_cor
#>    <int>   <dbl>        <dbl>         <dbl>              <dbl>    <dbl>
#> 1      1      24       0.0417       -0.0153            -0.0172  -0.0539
#> 2      2      24       0.0833       -0.0437            -0.0481   0.0179
#> 3      3      24       0.0417        0.0509             0.0387   0.0304
#> # i 6 more variables: era_diag_minus_off_same_block <dbl>,
#> #   era_diag_minus_off_diff_block <dbl>, era_lag_cor <dbl>,
#> #   geom_cor_partial <dbl>, geom_cor_run_partial <dbl>, geom_cor_xrun <dbl>
```

Key metrics include:

- `era_top1_acc` — top-1 encoding→retrieval accuracy for the item key.
- `era_diag_mean` — mean encoding–retrieval similarity on the diagonal.
- `era_diag_minus_off` — mean item-specific similarity: for each
  retrieval item, matched similarity minus its own mean similarity to
  nonmatching encoding items. With complete similarities this equals
  diagonal mean minus the grand off-diagonal mean. If pairwise
  similarities are missing, each item with an estimable contrast
  receives equal weight.
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

## 3. Relating item similarity to retrieval ratings

A retrieval rating can be related directly to an item-specific E–R score
computed in each region or searchlight sphere. By default, this is
matched similarity minus that retrieval item’s mean nonmatch similarity,
so every retrieval trial has its own local background. Missing ratings
do not remove items from `era_diag_mean`; each zero-order correlation
uses its own complete pairs. Set `era_association_score = "matched"`
only when the intended response is raw matched E–R similarity without
background correction.

``` r

set.seed(512)
toy$design$test_design$vividness <- sample(1:6, length(item_keys), replace = TRUE)
toy$design$test_design$vividness[c(4, 17)] <- NA_real_
toy$design$test_design$retrieval_run <- factor(
  rep(paste0("run_", 1:3), length.out = length(item_keys))
)
toy$design$test_design$trial_order <- seq_along(item_keys)
```

Use `era_correlates` for the unadjusted relationship. Use
`era_association` for the adjusted model, `era_effects` to select
directional one-degree-of-freedom effects, and `era_effects_block` for
omnibus tests of one or more terms:

``` r

era_assoc_ms <- era_rsa_model(
  dataset = toy$dataset,
  design = toy$design,
  key_var = ~ item,
  pairing = "one_to_one",
  era_simfun = "spearman",
  era_min_voxels = 3,
  era_components = "item",
  era_association_score = "item_specificity",
  era_correlates = ~ vividness,
  era_cor_method = "spearman",
  era_association = ~ vividness + retrieval_run + trial_order,
  era_effects = ~ vividness,
  era_effects_block = list(acquisition = ~ retrieval_run + trial_order),
  era_min_complete = 4
)

era_assoc_res <- run_regional(era_assoc_ms, region_mask)
era_assoc_res$performance_table[
  ,
  c("roinum", "era_diag_mean", "era_vividness_cor", "era_vividness_n",
    "era_assoc_part_r_vividness", "era_assoc_dr2_acquisition",
    "era_assoc_F_acquisition", "era_assoc_df1_acquisition",
    "era_assoc_n", "era_assoc_df_resid")
]
#> # A tibble: 3 x 10
#>   roinum era_diag_mean era_vividness_cor era_vividness_n era_assoc_part_r_vivi~1
#>    <int>         <dbl>             <dbl>           <dbl>                   <dbl>
#> 1      1       -0.0191           -0.0473              22                 -0.0332
#> 2      2       -0.0541           -0.0346              22                 -0.111 
#> 3      3        0.0504           -0.276               22                 -0.211 
#> # i abbreviated name: 1: era_assoc_part_r_vividness
#> # i 5 more variables: era_assoc_dr2_acquisition <dbl>,
#> #   era_assoc_F_acquisition <dbl>, era_assoc_df1_acquisition <dbl>,
#> #   era_assoc_n <dbl>, era_assoc_df_resid <dbl>
```

`era_vividness_cor` is the zero-order Spearman correlation with item
specificity. `era_assoc_part_r_vividness` is the signed semi-partial
correlation after adjusting for retrieval run and trial order. It is
invariant to linear rescaling of vividness. Adjustment-only terms remain
in the model but do not create unsigned or reference-dependent maps. A
multi-level factor needs an explicit one-degree-of-freedom contrast
before it can be requested in `era_effects`. Blocks do not have that
restriction: here, the `acquisition` block jointly removes the
multi-column retrieval-run factor and trial order.
`era_assoc_dr2_acquisition` is the block’s incremental R-squared,
`era_assoc_F_acquisition` is its partial F statistic, and
`era_assoc_df1_acquisition` is the fitted rank difference between the
full and reduced models.

The association model uses one joint complete-case set across similarity
and all regression variables, reported as `era_assoc_n`; every block
comparison reuses exactly those rows. It does not emit raw beta or
per-term unsigned R-squared maps. Before conventional second-level
analysis, Fisher-transform subject-level correlation or part-r maps
where that model’s assumptions call for it.

## 4. When do you need ERA partitioning?

The ERA-RSA summary above gives you direct item matching and a raw
geometry correlation.
[`era_partition_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
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
  include_procrustes = FALSE # keep the toy example fast; see note below
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
  [`naive_xdec_model()`](https://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md).
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

[`era_partition_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
can also include block, run, category, temporal, and custom nuisance
regressors through item-level vectors such as `item_block_enc`,
`item_block_ret`, `item_run_enc`, `item_run_ret`, `item_time_enc`,
`item_time_ret`, `item_category`, plus `first_order_nuisance` and
`second_order_nuisance`. By default, available item-level nuisance
groups are added automatically. Use `auto_nuisance = FALSE` to use only
custom nuisance regressors, or a character vector such as
`auto_nuisance = c("run", "time")` to keep selected groups. If you have
enough paired items, set `include_procrustes = TRUE` to add a
leakage-free orthogonal alignment decoder.

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

## 5. Searchlight ERA-RSA

We can also run ERA-RSA in a searchlight mode to obtain whole-brain maps
of the same metrics.

``` r

set.seed(456)

era_sl <- run_searchlight(
  era_assoc_ms,
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
#>   -  era_diag_mean  (Type:  DenseNeuroVol ) 
#>   -  era_diag_minus_off  (Type:  DenseNeuroVol ) 
#>   -  era_diag_minus_off_same_block  (Type:  DenseNeuroVol ) 
#>   -  era_diag_minus_off_diff_block  (Type:  DenseNeuroVol ) 
#>   -  era_lag_cor  (Type:  DenseNeuroVol ) 
#>   -  era_vividness_cor  (Type:  DenseNeuroVol ) 
#>   -  era_vividness_n  (Type:  DenseNeuroVol ) 
#>   -  era_assoc_part_r_vividness  (Type:  DenseNeuroVol ) 
#>   -  era_assoc_dr2_acquisition  (Type:  DenseNeuroVol ) 
#>   -  era_assoc_F_acquisition  (Type:  DenseNeuroVol ) 
#>   -  era_assoc_df1_acquisition  (Type:  DenseNeuroVol ) 
#>   -  era_assoc_n  (Type:  DenseNeuroVol ) 
#>   -  era_assoc_df_resid  (Type:  DenseNeuroVol )
```

The `searchlight_result` contains:

- `metrics`: names of the output maps (e.g., `era_diag_minus_off`,
  `era_vividness_cor`, and `era_assoc_part_r_vividness`),
- `results`: a list of `NeuroVol` maps, one per metric.

``` r

era_sl$metrics
#>  [1] "n_items"                       "era_diag_mean"                
#>  [3] "era_diag_minus_off"            "era_diag_minus_off_same_block"
#>  [5] "era_diag_minus_off_diff_block" "era_lag_cor"                  
#>  [7] "era_vividness_cor"             "era_vividness_n"              
#>  [9] "era_assoc_part_r_vividness"    "era_assoc_dr2_acquisition"    
#> [11] "era_assoc_F_acquisition"       "era_assoc_df1_acquisition"    
#> [13] "era_assoc_n"                   "era_assoc_df_resid"
```

We can save the searchlight maps using
[`save_results()`](https://bbuchsbaum.github.io/rMVPA/reference/save_results.md):

``` r

out_dir <- tempfile("era_rsa_sl_")
dir.create(out_dir, showWarnings = FALSE)

save_results(era_sl, out_dir, level = "standard")
list.files(file.path(out_dir, "maps"))
```

This will create one NIfTI file per requested metric (e.g.,
`era_diag_minus_off.nii.gz`, `era_vividness_cor.nii.gz`) that can be
viewed in your favorite neuroimaging software.

For association-focused HPC analyses, combine the item-only model above
with the shard backend:

``` r

era_sl_hpc <- run_searchlight(
  era_assoc_ms,
  radius = 2,
  method = "standard",
  engine = "era_rsa_fast",
  backend = "shard"
)
```

`era_components = "item"` skips top-1 identification and both
within-phase RDMs. For a one-to-one volumetric standard searchlight,
`engine = "auto"` (the default) selects the same dedicated ERA-RSA
engine shown explicitly above. It slices train/test matrices directly,
applies the established train-only feature filter and
center-preservation rule, and then calls the same
[`fit_roi.era_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/fit_roi.era_rsa_model.md)
computation used by the general-purpose iterator. With
`backend = "shard"`, workers slice shared matrices directly rather than
constructing and serializing an ROI object for every sphere. Request
`era_components = c("item", "identification")` when top-1 accuracy is
also needed, or retain the default components when geometry maps are
required; all three component families use the same dedicated extraction
engine. Repeated- item averaging, non-volumetric data, and
randomized/resampled searchlights continue through the general-purpose
iterator. Its existing compatibility key, `engine = "legacy"`, remains
available for an explicit reference run.

## 6. Adding Confounds and Lag Information

ERA-RSA can optionally incorporate item-level confounds and lag
variables. These are all defined at the **item** level, not the trial
level, and must align with the levels of `key_var`.

### Item-level confounds (`confound_rdms`)

`confound_rdms` is a named list of K×K matrices or `"dist"` objects
describing item-by-item nuisance structure (e.g., block/run/time), where
rows and columns correspond to item keys.

Run confounds can now be supplied in the simpler item-level form used by
the cross-run diagnostic. You pass `item_run_enc` and `item_run_ret`;
ERA-RSA derives the same-run RDMs internally for `geom_cor_run_partial`.

``` r

items <- levels(toy$design$train_design$item)

Mode <- function(x) { ux <- unique(x); ux[which.max(tabulate(match(x, ux)))] }
item_run_enc <- sapply(items, function(it) {
  Mode(toy$design$train_design$block_var[toy$design$train_design$item == it])
})
item_run_ret <- sapply(items, function(it) {
  Mode(toy$design$test_design$block_var[toy$design$test_design$item == it])
})
names(item_run_enc) <- names(item_run_ret) <- items

item_run_enc <- factor(paste0("enc_", item_run_enc))
item_run_ret <- factor(paste0("ret_", item_run_ret))
names(item_run_enc) <- names(item_run_ret) <- items
```

You can add temporal or block RDMs to `confound_rdms` and choose which
groups enter the general partial-geometry metric with `partial_against`.

``` r

time_enc <- setNames(seq_along(items), items)
time_enc_rdm <- abs(outer(time_enc, time_enc, "-"))
rownames(time_enc_rdm) <- colnames(time_enc_rdm) <- items

era_ms_conf <- era_rsa_model(
  dataset = toy$dataset,
  design  = toy$design,
  key_var = ~ item,
  phase_var = ~ block_var,
  confound_rdms = list(time_enc = time_enc_rdm),
  partial_against = c("run", "time"),
  item_run_enc = item_run_enc,
  item_run_ret = item_run_ret
)

era_conf_res <- run_regional(era_ms_conf, region_mask)
```

``` r

era_conf_res$performance_table[
  ,
  c("roinum", "geom_cor", "geom_cor_partial",
    "geom_cor_run_partial", "geom_cor_xrun", "beta_time_enc")
]
#> # A tibble: 3 x 6
#>   roinum geom_cor geom_cor_partial geom_cor_run_partial geom_cor_xrun
#>    <int>    <dbl>            <dbl>                <dbl>         <dbl>
#> 1      1 -0.0631          -0.0638              -0.0628             NA
#> 2      2  0.00263          0.00133              0.00252            NA
#> 3      3  0.0199           0.0159               0.0200             NA
#> # i 1 more variable: beta_time_enc <dbl>
```

The three geometry metrics have different interpretations:

- `geom_cor`: raw correlation between encoding and retrieval RDMs.
- `geom_cor_partial`: residualized geometry correlation after regressing
  out the nuisance groups selected by `partial_against`.
- `geom_cor_run_partial`: legacy run-only partial geometry correlation.
  It is still reported for continuity and is derived from `item_run_enc`
  / `item_run_ret` when explicit `confound_rdms$run_enc` and
  `confound_rdms$run_ret` are absent.
- `geom_cor_xrun`: raw geometry correlation restricted to item pairs
  that differ in both encoding and retrieval run.

For temporal confounds, a common choice is
`partial_against = c("run", "time")`. If you want only a temporal
residualization, use `partial_against = "time"` or the exact RDM name,
such as `partial_against = "time_enc"`.

### Whole-mask global nuisance

If you want to remove broad similarity structure shared across the whole
analysis mask, set `global_nuisance = TRUE`. ERA-RSA computes item-level
whole-mask RDMs once and adds them as `global_enc` and `global_ret`
nuisance RDMs. When `partial_against` is left at its default, enabling
`global_nuisance` makes the effective partial model
`c("run", "global")`; set `partial_against` explicitly if you want a
different nuisance set.

``` r

item_block <- factor(item_run_enc, levels = sort(unique(item_run_enc)))

era_ms_global <- era_rsa_model(
  dataset = toy$dataset,
  design  = toy$design,
  key_var = ~ item,
  phase_var = ~ block_var,
  item_block = item_block,
  item_run_enc = item_run_enc,
  item_run_ret = item_run_ret,
  global_nuisance = TRUE
)

era_global_res <- run_regional(era_ms_global, region_mask)
```

``` r

era_global_res$performance_table[
  ,
  c("roinum", "geom_cor_partial", "geom_cor_run_partial",
    "beta_global_enc", "beta_global_ret")
]
#> # A tibble: 3 x 5
#>   roinum geom_cor_partial geom_cor_run_partial beta_global_enc beta_global_ret
#>    <int>            <dbl>                <dbl>           <dbl>           <dbl>
#> 1      1           0.0152             -0.0628          -0.0989           1.16 
#> 2      2           0.0134              0.00252          0.0701           0.860
#> 3      3          -0.0260              0.0200           0.0711           0.850
```

### Block structure (`item_block`)

`item_block` encodes a per-item block/condition, typically derived from
`design$train_design$block_var`:

``` r

era_ms_block <- era_rsa_model(
  dataset    = toy$dataset,
  design     = toy$design,
  key_var    = ~ item,
  phase_var  = ~ block_var,
  item_block = item_block
)
```

This enables block-limited ERA contrasts:

- `era_diag_minus_off_same_block`: each retrieval item’s matched ERA
  minus its mean similarity to other encoding items in the same block,
  averaged equally over retrieval items.
- `era_diag_minus_off_diff_block`: each retrieval item’s matched ERA
  minus its mean similarity to encoding items in different blocks,
  averaged equally over retrieval items.

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
[`era_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md).

### Run and temporal nuisance in `era_partition_model()`

The variance-partition model uses the same item-level metadata but
separates first-order cross-state similarity from second-order geometry
preservation. When `auto_nuisance` includes `"run"`, the model adds
same-run nuisance regressors to both levels of the partition.

``` r

partition_nuis_ms <- era_partition_model(
  dataset = toy$dataset,
  design  = toy$design,
  key_var = ~ item,
  item_run_enc = item_run_enc,
  item_run_ret = item_run_ret,
  item_time_enc = time_enc,
  item_time_ret = time_enc + length(items),
  auto_nuisance = c("run", "time"),
  include_procrustes = FALSE
)

partition_nuis_res <- run_regional(partition_nuis_ms, region_mask)
```

``` r

partition_nuis_res$performance_table[
  ,
  c("roinum", "first_order_delta_r2", "second_order_delta_r2",
    "nuisance_first_order_n", "nuisance_second_order_n")
]
#> # A tibble: 3 x 5
#>   roinum first_order_delta_r2 second_order_delta_r2 nuisance_first_order_n
#>    <int>                <dbl>                 <dbl>                  <dbl>
#> 1      1             0.000244            0.00406                         4
#> 2      2             0.00198             0.00000176                      4
#> 3      3             0.00123             0.000252                        4
#> # i 1 more variable: nuisance_second_order_n <dbl>
```

Here the first-order nuisance set contains `same_run_cross`, `enc_time`,
`ret_time`, and `abs_lag`. The second-order nuisance set contains
`same_run_enc`, `same_run_ret`, `temporal_distance_enc`, and
`temporal_distance_ret`. Set `auto_nuisance = FALSE` when you want to
replace these defaults with custom kernels through
`first_order_nuisance` and `second_order_nuisance`.

Whole-mask nuisance is also available in the partition model. When
`auto_nuisance` includes `"global"`, the first-order partition receives
`global_cross`; the second-order partition receives `global_enc` and
`global_ret`.

``` r

partition_global_ms <- era_partition_model(
  dataset = toy$dataset,
  design  = toy$design,
  key_var = ~ item,
  global_nuisance = TRUE,
  auto_nuisance = "global",
  include_procrustes = FALSE
)

partition_global_res <- run_regional(partition_global_ms, region_mask)
```

``` r

partition_global_res$performance_table[
  ,
  c("roinum", "first_order_delta_r2", "second_order_delta_r2",
    "nuisance_first_order_n", "nuisance_second_order_n")
]
#> # A tibble: 3 x 5
#>   roinum first_order_delta_r2 second_order_delta_r2 nuisance_first_order_n
#>    <int>                <dbl>                 <dbl>                  <dbl>
#> 1      1            0.0000202              0.000128                      1
#> 2      2            0.00113                0.000139                      1
#> 3      3            0.00236                0.000493                      1
#> # i 1 more variable: nuisance_second_order_n <dbl>
```

## 7. Summary

- [`era_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md)
  provides a unified framework for:
  - cross-decoding between encoding and retrieval, and
  - comparing encoding and retrieval representational geometries, and
  - mapping zero-order and adjusted directional associations between
    item-specific matched-minus-nonmatch similarity and retrieval
    variables.
- [`era_partition_model()`](https://bbuchsbaum.github.io/rMVPA/reference/era_partition_model.md)
  separates direct item transfer, variance uniquely explained by
  same-item similarity, and variance explained by preserved second-order
  geometry.
- It integrates naturally with:
  - [`run_regional()`](https://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
    for ROI-based analyses, and
  - [`run_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
    for whole-brain mapping.
- Outputs are standard rMVPA result objects, so the same tooling
  ([`save_results()`](https://bbuchsbaum.github.io/rMVPA/reference/save_results.md),
  plotting, etc.) applies as for other models such as
  [`naive_xdec_model()`](https://bbuchsbaum.github.io/rMVPA/reference/naive_xdec_model.md)
  and
  [`remap_rrr_model()`](https://bbuchsbaum.github.io/rMVPA/reference/remap_rrr_model.md).

ERA-RSA is particularly useful when you want to go beyond simple
accuracy and ask **how similar the geometry of neural representations is
across phases**, and how that similarity varies across brain regions.
