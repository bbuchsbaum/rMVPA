# Feature-Based RSA

## What problem does Feature-Based RSA solve?

You already have a feature space for your stimuli — semantic embeddings,
CNN layer activations, behavioural ratings, model logits. The question
is not “do two RDMs correlate?” but **“does this brain region encode
these features well enough that I can predict held-out neural patterns
from them?”**

Feature-Based RSA fits a cross-validated regression from the feature
matrix `F` (trials × features) to the neural pattern matrix `X` (trials
× voxels), evaluates on held-out trials, and reports several
reconstruction metrics. It sits between standard RSA (correlate two
RDMs) and full encoding/decoding (fit voxelwise GLMs): you get an
RDM-level summary *and* trial-level pattern fidelity, with
regularisation and component selection built in.

## When to use it

Use Feature-Based RSA when:

- You have a **named, fixed feature space** (rows = trials, columns =
  features) and want to ask whether an ROI encodes those specific
  dimensions.
- You care about **per-trial pattern reconstruction**, not just whether
  two RDMs correlate — for example, you want to test whether semantic
  embeddings predict the *exact* held-out trial pattern, not just
  average geometry.
- You want **regularisation** (PLS, PCR, or elastic net) because your
  feature space is high-dimensional, correlated, or has more features
  than trials.
- You want **multiple performance metrics** out of one fit — pattern
  correlation, pattern discrimination, RDM correlation, voxel
  correlation, R² — rather than choosing one a priori.

If you only have a model RDM (no feature matrix) and want a per-trial
similarity score with built-in across-block masking, see
[`vignette("Vector_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Vector_RSA.md).
If you want to test multiple model RDMs as regressors and read out
coefficients, see
[`vignette("RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/RSA.md).

## Inputs and outputs

| Component | What it is |
|----|----|
| `F` (or `S`) | Feature matrix `n_trials × n_features` — *or* a similarity matrix `S` that gets eigen-decomposed into a feature basis. |
| `labels` | Vector of length `n_trials` naming each row of `F`. |
| `mvpa_dataset` | Neural data `X`: `n_trials × n_voxels` per ROI/searchlight. |
| `crossval` | A cross-validation spec. Required (or pass `block_var` to the design and let it build a blocked CV). |
| Output | Per-ROI metrics including `pattern_correlation`, `pattern_discrimination`, `pattern_rank_percentile`, `rdm_correlation`, `voxel_correlation`, `mse`, `r_squared`, `mean_voxelwise_temporal_cor`, and `ncomp` (components actually used). |

## Minimal example

We will:

1.  Generate a small synthetic dataset where the brain patterns *are*
    driven by a known feature matrix.
2.  Build a `feature_rsa_design` from that feature matrix.
3.  Fit
    [`feature_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
    with PLS regression and inspect the per-region performance.

``` r

set.seed(2026)

# Synthetic neural data: 50 trials in 4 blocks, small mask.
sim <- gen_sample_dataset(D = c(6, 6, 6), nobs = 50, blocks = 4, nlevels = 2)

# Feature matrix shared across the simulation: 5 features driven by 2 latent factors.
n_trials   <- 50
n_features <- 5
Z <- matrix(rnorm(n_trials * 2), n_trials, 2)
B <- matrix(rnorm(2 * n_features), 2, n_features)
F <- 0.7 * (Z %*% B) + 0.3 * matrix(rnorm(n_trials * n_features), n_trials, n_features)
F <- base::scale(F)
colnames(F) <- paste0("feat_", seq_len(n_features))
labels <- paste0("trial_", seq_len(n_trials))
```

Inject a weak signal driven by the *same* latent factors `Z` into the
brain data, so that a feature-RSA model has something to recover:

``` r

mask_vol  <- sim$dataset$mask
mask_idx  <- which(as.logical(neuroim2::values(mask_vol)))

datamat <- do.call(cbind, lapply(
  neuroim2::vols(sim$dataset$train_data),
  function(v) as.numeric(v[mask_idx])
))

p1 <- rnorm(length(mask_idx)); p1 <- p1 / sd(p1)
p2 <- rnorm(length(mask_idx)); p2 <- p2 / sd(p2)
datamat <- datamat + p1 %*% t(0.5 * Z[, 1]) + p2 %*% t(0.5 * Z[, 2])

train_vec <- neuroim2::SparseNeuroVec(
  datamat, neuroim2::space(sim$dataset$train_data),
  mask = as.logical(neuroim2::values(mask_vol))
)
dset <- mvpa_dataset(train_vec, mask = mask_vol)
```

Now build the design and the model:

``` r

design <- feature_rsa_design(
  F        = F,
  labels   = labels,
  max_comps = 2          # cap PLS/PCR components
)

cv <- blocked_cross_validation(sim$design$block_var)

mod <- feature_rsa_model(
  dataset  = dset,
  design   = design,
  method   = "pls",      # PLS regression: F -> X
  crossval = cv
)
```

A regional run that splits the mask into a few ROIs:

``` r

nvox <- sum(mask_vol)
region_mask <- neuroim2::NeuroVol(
  sample(1:3, size = nvox, replace = TRUE),
  neuroim2::space(mask_vol),
  indices = which(mask_vol > 0)
)

res <- run_regional(mod, region_mask)
res$performance_table
#> # A tibble: 3 × 10
#>   roinum pattern_correlation pattern_discrimination pattern_rank_percentile
#>    <int>               <dbl>                  <dbl>                   <dbl>
#> 1      1               0.385                  0.344                   0.788
#> 2      2               0.405                  0.374                   0.798
#> 3      3               0.361                  0.339                   0.774
#> # ℹ 6 more variables: rdm_correlation <dbl>, voxel_correlation <dbl>,
#> #   mse <dbl>, r_squared <dbl>, mean_voxelwise_temporal_cor <dbl>, ncomp <dbl>
```

Each row is one ROI. `pattern_correlation` and `rdm_correlation` are
usually the most diagnostic columns: they answer “did the model produce
a sensible held-out pattern?” and “did the predicted pattern have the
right pairwise geometry?” respectively. `ncomp` records how many PLS
components were actually used (controlled by `ncomp_selection`).

## Choosing a method

[`feature_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
supports three estimators. The right choice depends on the feature
space:

| Method | Best when | Notes |
|----|----|----|
| `"pls"` (default) | Features are correlated; you want supervised dimensionality reduction. | Uses [`pls::plsr`](https://khliland.github.io/pls/reference/mvr.html). Fits up to `max_comps` components. |
| `"pca"` | You want to reduce features unsupervised first, then regress on PCs. | Principal component regression via [`pls::pcr`](https://khliland.github.io/pls/reference/mvr.html). |
| `"glmnet"` | High-dimensional or sparse features; you want shrinkage. | Multivariate Gaussian elastic net via `glmnet`. Tune with `alpha`, `lambda`, `cv_glmnet`. |

For PLS and PCR, `ncomp_selection` decides how many of the fitted
components are used at prediction time:

- `"loo"` (default): leave-one-out CV, picks the fewest components
  within one SE of the minimum RMSEP.
- `"pve"`: stop when cumulative explained variance crosses
  `pve_threshold` (default 0.9).
- `"max"`: always use `max_comps` components.

``` r

methods <- c("pls", "pca", "glmnet")
summary_rows <- lapply(methods, function(m) {
  spec <- feature_rsa_model(dataset = dset, design = design, method = m, crossval = cv)
  out  <- run_regional(spec, region_mask)
  cbind(method = m, out$performance_table)
})
do.call(rbind, summary_rows)
#>   method roinum pattern_correlation pattern_discrimination
#> 1    pls      1           0.3846940              0.3441974
#> 2    pls      2           0.4045292              0.3739388
#> 3    pls      3           0.3612074              0.3388624
#> 4    pca      1           0.3805947              0.3398303
#> 5    pca      2           0.4050113              0.3740592
#> 6    pca      3           0.3576422              0.3349944
#> 7 glmnet      1           0.3719819              0.3342702
#> 8 glmnet      2           0.3859255              0.3572923
#> 9 glmnet      3           0.3374512              0.3165107
#>   pattern_rank_percentile rdm_correlation voxel_correlation      mse r_squared
#> 1               0.7877551       0.6901634         0.4523082 1.143468 0.1926938
#> 2               0.7975510       0.7482323         0.4646539 1.074230 0.2044189
#> 3               0.7738776       0.6653580         0.4313721 1.197495 0.1720092
#> 4               0.7848980       0.6827079         0.4468376 1.150273 0.1878893
#> 5               0.7991837       0.7515874         0.4648338 1.072501 0.2056996
#> 6               0.7734694       0.6618135         0.4279024 1.200183 0.1701507
#> 7               0.8016327       0.7153794         0.4375939 1.176173 0.1696034
#> 8               0.8020408       0.7590258         0.4450374 1.114066 0.1749165
#> 9               0.7751020       0.6942123         0.4033873 1.258269 0.1299875
#>   mean_voxelwise_temporal_cor ncomp
#> 1                   0.3148281     2
#> 2                   0.3338867     2
#> 3                   0.2961254     2
#> 4                   0.3096119     2
#> 5                   0.3341077     2
#> 6                   0.2929064     2
#> 7                   0.3121517     5
#> 8                   0.3287285     5
#> 9                   0.2806889     5
```

## How it differs from standard RSA

|  | Standard RSA (`rsa_model`) | Feature-Based RSA (`feature_rsa_model`) |
|----|----|----|
| Input | One or more **model RDMs** | A **feature matrix** `F` (or symmetric similarity `S`) |
| Comparison | Correlation/regression on **vectorised RDMs** | Cross-validated regression on **trial-level patterns** |
| Cross-validation | Optional (block exclusion) | Required (`crossval`) — predicts held-out trials |
| Output | One coefficient or correlation per RDM | Multiple metrics per ROI: pattern, RDM, voxel, R² |
| Good for | Testing competing geometric hypotheses | Asking *“does this region encode these features?”* |

Standard RSA is lighter and more interpretable when your hypothesis is
about geometry. Feature-Based RSA is the right tool when you have an
actual feature space and want to know whether the brain tracks it well
enough to *predict* the next held-out trial.

## What’s next

- [`vignette("Vector_RSA")`](http://bbuchsbaum.github.io/rMVPA/articles/Vector_RSA.md)
  — per-trial RSA scores with built-in across-block masking, when you
  have a single reference RDM.
- [`vignette("Feature_RSA_Advanced_Workflows")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Advanced_Workflows.md)
  — workflow extensions: returning predicted RDM vectors, cross-ROI
  representational connectivity from feature-RSA fits, multi-scenario
  evaluation.
- [`vignette("Feature_RSA_Connectivity")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Connectivity.md)
  and
  [`vignette("Feature_RSA_Domain_Adaptation")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Domain_Adaptation.md)
  — connectivity and cross-state extensions of the same model.
- [`?feature_rsa_design`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md),
  [`?feature_rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
  — full argument reference.
