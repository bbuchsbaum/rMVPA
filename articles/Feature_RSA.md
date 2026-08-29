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
[`vignette("Vector_RSA")`](https://bbuchsbaum.github.io/rMVPA/articles/Vector_RSA.md).
If you want to test multiple model RDMs as regressors and read out
coefficients, see
[`vignette("RSA")`](https://bbuchsbaum.github.io/rMVPA/articles/RSA.md).

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
    [`feature_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
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
  max_comps = 2,         # cap PLS/PCR components
  block_var = sim$design$block_var
)

cv <- blocked_cross_validation(sim$design$block_var)

mod <- feature_rsa_model(
  dataset  = dset,
  design   = design,
  method   = "pls",      # PLS regression: F -> X
  ncomp_selection = "blocked",
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
#> 1      1               0.285                  0.228                   0.669
#> 2      2               0.325                  0.280                   0.7  
#> 3      3               0.251                  0.216                   0.654
#> # ℹ 6 more variables: rdm_correlation <dbl>, voxel_correlation <dbl>,
#> #   mse <dbl>, r_squared <dbl>, mean_voxelwise_temporal_cor <dbl>, ncomp <dbl>
```

Each row is one ROI. `pattern_correlation` and `rdm_correlation` are
usually the most diagnostic columns: they answer “did the model produce
a sensible held-out pattern?” and “did the predicted pattern have the
right pairwise geometry?” respectively. `ncomp` records how many PLS
components were actually used (controlled by `ncomp_selection`).

## Choosing a method

[`feature_rsa_model()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
supports three estimators. The right choice depends on the feature
space:

| Method | Best when | Notes |
|----|----|----|
| `"pls"` (default) | Features are correlated; you want supervised dimensionality reduction. | Uses the configured `pls` PLS numerical kernel. Fits up to `max_comps` components. |
| `"pca"` | You want to reduce features unsupervised first, then regress on PCs. | Uses the configured `pls` PCR numerical kernel (SVD-PCR by default). |
| `"glmnet"` | High-dimensional or sparse features; you want shrinkage. | Multivariate Gaussian elastic net via `glmnet`. Tune with `alpha`, `lambda`, `cv_glmnet`. |

## How should you select the number of components?

PLS and PCR fit up to `max_comps`; `ncomp_selection` decides how many
are used for each outer-fold prediction. This tuning happens **inside**
the outer cross-validation loop, so held-out outer observations never
choose their own component count.

| Selection | Inner work per outer training fold | What it estimates | When to use it |
|----|---:|----|----|
| `"blocked"` | One fit per training block, plus the final fit | Leave-one-block-out prediction error; each block contributes one MSE | Preferred when runs, sessions, subjects, or another genuine independence unit are available |
| `"loo"` (default) | One fit per training observation, plus the final fit | Leave-one-observation-out prediction error | Backward-compatible choice for exchangeable, independent rows; expensive for long designs |
| `"pve"` | One fit | Cumulative predictor variance explained, stopping at `pve_threshold` (default 0.9) | A fast heuristic when held-out tuning is not required |
| `"max"` | One fit | No tuning; always use the available `max_comps` | A fixed-complexity analysis chosen in advance |

The `"blocked"` and `"loo"` selectors use the same one-standard-error
rule: choose the fewest components whose mean segment MSE is within one
standard error of the minimum. They differ in the unit treated as
exchangeable. If observations within a run are autocorrelated, LOO can
be both unnecessarily slow and scientifically optimistic; blocking on
run makes the validation unit match the acquisition structure. For
`"blocked"`, feature and neural centering and scaling are re-estimated
from each inner training split, so the held-out block does not set
preprocessing parameters. With unequal block sizes, blocks receive equal
weight rather than observations receiving equal weight.

The default remains `"loo"` for compatibility. On blocked neuroimaging
data, pass the same observation-aligned `block_var` to
[`feature_rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md)
and use `ncomp_selection = "blocked"`, as in the minimal example. Every
outer training fold must retain at least two non-empty blocks; rMVPA
fails the fold rather than silently reverting to a different selector.

The implementation streams validation errors and retains only the
compact coefficient path needed for prediction. This reduces transient
memory but does not change the number of fits implied by the selector.
In particular, exact LOO remains intrinsically proportional to the
number of training observations. Parallel workers keep one target matrix
plus fold indices instead of serializing overlapping target slices for
every fold.

``` r

methods <- c("pls", "pca", "glmnet")
summary_rows <- lapply(methods, function(m) {
  spec <- feature_rsa_model(
    dataset = dset, design = design, method = m, crossval = cv,
    ncomp_selection = "blocked"
  )
  out  <- run_regional(spec, region_mask)
  cbind(method = m, out$performance_table)
})
do.call(rbind, summary_rows)
#>   method roinum pattern_correlation pattern_discrimination
#> 1    pls      1           0.2850026              0.2275442
#> 2    pls      2           0.3247986              0.2802777
#> 3    pls      3           0.2512543              0.2157884
#> 4    pca      1           0.2835266              0.2262611
#> 5    pca      2           0.3254417              0.2808132
#> 6    pca      3           0.2477400              0.2120167
#> 7 glmnet      1           0.3719819              0.3342702
#> 8 glmnet      2           0.3859255              0.3572923
#> 9 glmnet      3           0.3374512              0.3165107
#>   pattern_rank_percentile rdm_correlation voxel_correlation      mse r_squared
#> 1               0.6689796       0.5496127         0.3649092 1.248759 0.1183568
#> 2               0.7000000       0.6694226         0.4053769 1.142053 0.1541890
#> 3               0.6542857       0.5335336         0.3475725 1.298152 0.1024111
#> 4               0.6685714       0.5494108         0.3604843 1.254104 0.1145828
#> 5               0.7024490       0.6688983         0.4043075 1.143148 0.1533781
#> 6               0.6538776       0.5240422         0.3447385 1.300160 0.1010227
#> 7               0.8016327       0.7153794         0.4375939 1.176173 0.1696034
#> 8               0.8020408       0.7590258         0.4450374 1.114066 0.1749165
#> 9               0.7751020       0.6942123         0.4033873 1.258269 0.1299875
#>   mean_voxelwise_temporal_cor ncomp
#> 1                   0.2293049     2
#> 2                   0.2684918     2
#> 3                   0.2130810     2
#> 4                   0.2240035     2
#> 5                   0.2669684     2
#> 6                   0.2096067     2
#> 7                   0.3121517     5
#> 8                   0.3287285     5
#> 9                   0.2806889     5
```

For a reproducible performance characterization, run
`inst/benchmarks/feature_rsa_hotpaths.R` from a source checkout. It
compares the compact kernels with the former formula-based
implementation, validates predictions and segment errors before timing,
and records absolute repeated times rather than enforcing a
machine-specific speed claim in ordinary tests.

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

- [`vignette("Vector_RSA")`](https://bbuchsbaum.github.io/rMVPA/articles/Vector_RSA.md)
  — per-trial RSA scores with built-in across-block masking, when you
  have a single reference RDM.
- [`vignette("Feature_RSA_Advanced_Workflows")`](https://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Advanced_Workflows.md)
  — workflow extensions: returning predicted RDM vectors, cross-ROI
  representational connectivity from feature-RSA fits, multi-scenario
  evaluation.
- [`vignette("Feature_RSA_Connectivity")`](https://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Connectivity.md)
  and
  [`vignette("Feature_RSA_Domain_Adaptation")`](https://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Domain_Adaptation.md)
  — connectivity and cross-state extensions of the same model.
- [`?feature_rsa_design`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md),
  [`?feature_rsa_model`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
  — full argument reference.
