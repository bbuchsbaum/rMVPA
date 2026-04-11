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
  D            = c(6, 6, 6),  # modest volume for speed
  nobs         = 48,
  nlevels      = 3,           # 3 item categories
  blocks       = 3,           # 3 runs
  external_test = TRUE        # required for encoding → retrieval
)

# Add an item key column to both train and test designs
toy$design$train_design$item <- toy$design$train_design$Y
toy$design$test_design$item  <- toy$design$test_design$Ytest

str(toy$design$train_design[item = 1:6, ])
```

    tibble [48 x 4] (S3: tbl_df/tbl/data.frame)
     $ Y        : Factor w/ 3 levels "a","b","c": 1 3 1 3 3 1 2 3 1 3 ...
     $ block_var: int [1:48] 1 1 1 1 1 1 1 1 1 1 ...
     $ .rownum  : int [1:48] 1 2 3 4 5 6 7 8 9 10 ...
     $ item     : Factor w/ 3 levels "a","b","c": 1 3 1 3 3 1 2 3 1 3 ...

Each row in `train_design` and `test_design` corresponds to a trial. The
`item` column links encoding and retrieval presentations of the same
stimulus.

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
  phase_var = ~ block_var,     # phase label (not critical in external-test path)
  distfun   = cordist("pearson"),
  rsa_simfun = "spearman"
)

era_ms
```

     Model Specification 

    - Summary 
      - Name:  era_rsa_model 
      - Primary Class:  era_rsa_model 
      - Class Chain:  era_rsa_model, model_spec, list 
      - Return Predictions:  FALSE 
      - Compute Performance:  TRUE 
      - Has Test Set:  TRUE 
      - Basis Count:  1 
    - Parameters 
      - key_var: ~item
      - phase_var: ~block_var
      - encoding_level: 1
      - retrieval_level: 2
      - distfun: cordist / pearson
      - rsa_simfun: spearman
      - include_diag: TRUE

     MVPA Dataset 

    - Training Data 
      - Dimensions:  6 x 6 x 6 x 48 observations 
      - Type:  DenseNeuroVec 
    - Test Data 
      - Dimensions:  6 x 6 x 6 x 48 observations 
      - Type:  DenseNeuroVec 
    - Mask Information 
      - Areas:  TRUE : 216 
      - Active voxels/vertices:  216 


     MVPA Design 

    - Training Data 
      - Observations:  48 
      - Response Type:  Factor
      - Levels:  a, b, c 
      - Class Distribution:  a: 16, b: 16, c: 16 
    - Test Data 
      - Observations:  48 
      - Class Distribution:  a: 16, b: 16, c: 16 
    - Structure 
      - Blocking:  Present
      - Number of Blocks:  3 
      - Mean Block Size:  16  (SD:  0 ) 
      - Split Groups:  None 

``` r
era_res <- run_regional(era_ms, region_mask)
era_res
```

     Regional Analysis Results 

    - Summary 
      - Model:  era_rsa_model 
      - Regions with Results:  3 
      - Metrics:  n_items, era_top1_acc, era_diag_mean, era_diag_minus_off, geom_cor, era_diag_minus_off_same_block ... 
      - Output Maps:  n_items, era_top1_acc, era_diag_mean, era_diag_minus_off, geom_cor, era_diag_minus_off_same_block ... 

The `performance_table` contains one row per region and several ERA-RSA
metrics.

``` r
head(era_res$performance_table)
```

    # A tibble: 3 x 11
      roinum n_items era_top1_acc era_diag_mean era_diag_minus_off geom_cor
       <int>   <dbl>        <dbl>         <dbl>              <dbl>    <dbl>
    1      1       3        0.333       0.0555              0.0748     -1  
    2      2       3        0.667       0.00944             0.0428      0.5
    3      3       3        0.333       0.0473              0.0813     -0.5
    # i 5 more variables: era_diag_minus_off_same_block <dbl>,
    #   era_diag_minus_off_diff_block <dbl>, era_lag_cor <dbl>,
    #   geom_cor_run_partial <dbl>, geom_cor_xrun <dbl>

Key metrics include:

- `era_top1_acc` — top-1 encoding→retrieval accuracy for the item key.
- `era_diag_mean` — mean encoding–retrieval similarity on the diagonal.
- `era_diag_minus_off` — diagonal minus off-diagonal similarity.
- `geom_cor` — correlation between encoding and retrieval RDMs.

These quantify both cross-decoding performance and representational
geometry alignment between phases.

## 3. Searchlight ERA-RSA

We can also run ERA-RSA in a searchlight mode to obtain whole-brain maps
of the same metrics.

``` r
set.seed(456)

era_sl <- run_searchlight(
  era_ms,
  radius = 3,         # searchlight radius in voxels
  method = "standard" # ERA-RSA currently uses standard searchlight
)

era_sl
```

     Searchlight Analysis Results 

    - Coverage 
      - Voxels/Vertices in Mask:  216 
      - Voxels/Vertices with Results:  216 
    - Output Maps (Metrics) 
      -  n_items  (Type:  DenseNeuroVol ) 
      -  era_top1_acc  (Type:  DenseNeuroVol ) 
      -  era_diag_mean  (Type:  DenseNeuroVol ) 
      -  era_diag_minus_off  (Type:  DenseNeuroVol ) 
      -  geom_cor  (Type:  DenseNeuroVol ) 
      -  era_diag_minus_off_same_block  (Type:  DenseNeuroVol ) 
      -  era_diag_minus_off_diff_block  (Type:  DenseNeuroVol ) 
      -  era_lag_cor  (Type:  DenseNeuroVol ) 
      -  geom_cor_run_partial  (Type:  DenseNeuroVol ) 
      -  geom_cor_xrun  (Type:  DenseNeuroVol ) 

The `searchlight_result` contains:

- `metrics`: names of the output maps (e.g., `geom_cor`,
  `era_top1_acc`),
- `results`: a list of `NeuroVol` maps, one per metric.

``` r
era_sl$metrics
```

     [1] "n_items"                       "era_top1_acc"                 
     [3] "era_diag_mean"                 "era_diag_minus_off"           
     [5] "geom_cor"                      "era_diag_minus_off_same_block"
     [7] "era_diag_minus_off_diff_block" "era_lag_cor"                  
     [9] "geom_cor_run_partial"          "geom_cor_xrun"                

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

## 4. Adding Confounds and Lag Information

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

## 5. Summary

- [`era_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/era_rsa_model.md)
  provides a unified framework for:
  - cross-decoding between encoding and retrieval, and
  - comparing encoding and retrieval representational geometries.
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
