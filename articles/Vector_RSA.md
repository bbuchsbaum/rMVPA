# Vector-Based RSA

## What problem does Vector-Based RSA solve?

Standard RSA correlates the lower triangle of a neural RDM with the
lower triangle of a model RDM and gives you **one number per ROI**.
Vector-Based RSA does something subtly different: for **each trial `i`**
it correlates that trial’s *row* of the neural distance matrix with the
same row of a reference distance matrix — over the trials in *other*
blocks only — and returns one similarity score per trial. The mean
across trials becomes the ROI’s `rsa_score`.

This per-trial framing buys you three things:

1.  **Built-in across-block masking.** Within-run pairs (which suffer
    from temporal autocorrelation, scanner drift, breathing, etc.) are
    excluded automatically: for trial `i`, only columns where
    `block_var != block_var[i]` are used.
2.  **Per-observation scores.** With `return_predictions = TRUE` you get
    the full vector of trial-level RSA values, useful for trial-level
    regressions, item-level analyses, or single-trial QC.
3.  **A clean permutation null.** Set `nperm > 0` and a row-permutation
    null is computed *inside the ROI*, returning `p_rsa_score` and
    `z_rsa_score` alongside the mean.

## When to use it

Use Vector-Based RSA when:

- You have **one reference RDM** (from features, behaviour, a model, or
  a prior dataset) and want a single, interpretable RSA score per
  ROI/searchlight centre.
- Your design has **runs/blocks** and you want within-block comparisons
  excluded — this is the default, not an opt-in.
- You want **per-trial RSA scores** so you can correlate them with
  behaviour or a single-trial covariate.
- You want a **permutation-based p-value and z-score** without writing
  the permutation loop yourself.

If you instead have a feature matrix and want regression-style metrics
(pattern correlation, R²) under cross-validation, see
[`vignette("Feature_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA.md).
If you want to fit multiple model RDMs as competing predictors, see
[`vignette("RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/RSA.md).

## Inputs and outputs

| Component | What it is |
|----|----|
| `D` | Reference distance matrix, **with `rownames(D)` set**. Rows/columns are the *unique* labels (one row per condition or stimulus). Need not be the same length as the number of trials. |
| `labels` | Length `n_trials` vector. Each entry must appear in `rownames(D)`. The design expands `D` to an `n_trials × n_trials` matrix (`Dexpanded`) by row matching. |
| `block_var` | Length `n_trials` vector of block IDs (e.g. run number). Trials sharing a block are excluded from each other’s correlation neighbourhood. |
| `mvpa_dataset` | Neural data: `n_trials × n_voxels` per ROI. |
| `distfun` | A distance object from [`create_dist()`](http://bbuchsbaum.github.io/rMVPA/reference/create_dist.md) (e.g. [`cordist()`](http://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md), [`eucdist()`](http://bbuchsbaum.github.io/rMVPA/reference/distance-constructors.md)) used on the *neural* data. |
| `rsa_simfun` | `"pearson"` (default) or `"spearman"` — the correlation between the neural and reference distance vectors. |
| Output | Per-ROI `rsa_score` (mean across trials). With `nperm > 0`: also `p_rsa_score`, `z_rsa_score`. With `return_predictions = TRUE`: a `prediction_table` containing the per-trial scores. |

## Minimal example

``` r

set.seed(2026)

# Synthetic neural data: 50 trials in 5 blocks, small mask.
sim <- gen_sample_dataset(D = c(6, 6, 6), nobs = 50, blocks = 5, nlevels = 2)

# A reference distance matrix derived from a feature space.
n_trials   <- 50
n_features <- 5
F          <- matrix(rnorm(n_trials * n_features), n_trials, n_features)
labels     <- paste0("trial_", seq_len(n_trials))
rownames(F) <- labels

D <- as.matrix(dist(F))
rownames(D) <- colnames(D) <- labels
```

Build the design and the model:

``` r

design <- vector_rsa_design(
  D         = D,
  labels    = labels,
  block_var = sim$design$block_var
)

mod <- vector_rsa_model(
  dataset    = sim$dataset,
  design     = design,
  distfun    = cordist(),  # neural distance: 1 - Pearson correlation
  rsa_simfun = "pearson"
)
```

Run regionally over a few ROIs:

``` r

mask_vol <- sim$dataset$mask
nvox     <- sum(mask_vol)
region_mask <- neuroim2::NeuroVol(
  sample(1:3, size = nvox, replace = TRUE),
  neuroim2::space(mask_vol),
  indices = which(mask_vol > 0)
)

res <- run_regional(mod, region_mask)
res$performance_table
#> # A tibble: 3 × 2
#>   roinum rsa_score
#>    <int>     <dbl>
#> 1      1 -0.00647 
#> 2      2 -0.0121  
#> 3      3 -0.000355
```

Each row is one ROI; `rsa_score` is the average across-block
second-order correlation.

## What each row actually computes

For trial `i`:

    neural_dist_row_i  =  pairwise_dist(distfun, X)[i, valid]
    reference_row_i    =  Dexpanded[i, valid]
    score_i            =  cor(neural_dist_row_i, reference_row_i, method = rsa_simfun)

where `valid = which(block_var != block_var[i])`. The ROI’s `rsa_score`
is `mean(score_i, na.rm = TRUE)`. If a trial has no valid out-of-block
partners, its score is `NA` and it contributes nothing to the mean.

## Per-trial scores

If you want the full vector of scores back from `run_regional` or
`run_searchlight`, opt in:

``` r

mod_pred <- vector_rsa_model(
  dataset            = sim$dataset,
  design             = design,
  distfun            = cordist(),
  rsa_simfun         = "pearson",
  return_predictions = TRUE
)

res_pred <- run_regional(mod_pred, region_mask)
head(res_pred$prediction_table)
#> # A tibble: 6 × 3
#>   roinum rsa_score observation_index
#>    <int>     <dbl>             <int>
#> 1      1    0.0886                 1
#> 2      1   -0.0958                 2
#> 3      1   -0.384                  3
#> 4      1   -0.170                  4
#> 5      1    0.137                  5
#> 6      1   -0.187                  6
```

The `prediction_table` is one row per trial-per-ROI, with the
trial-level `rsa_score`. This lets you correlate RSA strength against
trial-level behavioural covariates without re-running the analysis.

## Permutation testing

Setting `nperm > 0` permutes the rows of the ROI’s data matrix and
recomputes the mean second-order similarity under the null. The output
adds `p_rsa_score` and `z_rsa_score`:

``` r

mod_perm <- vector_rsa_model(
  dataset    = sim$dataset,
  design     = design,
  distfun    = cordist(),
  rsa_simfun = "pearson",
  nperm      = 50            # use a larger value (e.g. 1000+) in real analyses
)

res_perm <- run_regional(mod_perm, region_mask)
res_perm$performance_table
#> # A tibble: 3 × 4
#>   roinum rsa_score p_rsa_score z_rsa_score
#>    <int>     <dbl>       <dbl>       <dbl>
#> 1      1 -0.00647        0.510     -0.206 
#> 2      2 -0.0121         0.608     -0.334 
#> 3      3 -0.000355       0.510      0.0149
```

Permutation is done *within each ROI* — it shuffles which neural row
goes with which reference row, holding the block structure of the
reference fixed. With `save_distributions = TRUE` the full null
distributions are also returned.

## How it differs from standard RSA

|  | Standard RSA (`rsa_model`) | Vector-Based RSA (`vector_rsa_model`) |
|----|----|----|
| Granularity | One scalar per ROI from the **upper-triangle vector** of pair dissimilarities | One scalar **per trial**, averaged across trials per ROI |
| Block handling | Optional, via `keep_intra_run = FALSE` | **Built-in** — within-block columns excluded for every trial |
| Multiple model RDMs | Supported as a regression formula | Single reference RDM at a time |
| Output | Coefficients or correlations per model RDM | `rsa_score` + optional `p_rsa_score`, `z_rsa_score`, per-trial `prediction_table` |
| Best for | Comparing competing geometric hypotheses | Trial-level RSA, single-RDM tests with strong block control |

## What’s next

- [`vignette("Feature_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA.md)
  — when you have a feature matrix rather than a single reference RDM.
- [`vignette("RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/RSA.md)
  — the standard whole-RDM workflow with multiple model RDMs and
  `keep_intra_run`.
- [`?vector_rsa_design`](http://bbuchsbaum.github.io/rMVPA/reference/vector_rsa_design.md),
  [`?vector_rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/vector_rsa_model.md),
  [`?second_order_similarity`](http://bbuchsbaum.github.io/rMVPA/reference/second_order_similarity.md)
  — argument reference and the underlying scoring function.
