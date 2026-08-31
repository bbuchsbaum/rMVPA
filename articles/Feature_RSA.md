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
- You want **regularisation** (PLS, PCR, fast ridge, or elastic net)
  because your feature space is high-dimensional, correlated, or has
  more features than trials.
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
| Output | Per-ROI reconstruction and geometry metrics. PLS/PCR and glmnet also report `ncomp` (a historical nonzero-count proxy for glmnet); ridge reports its selected-penalty, effective-dimension, and grid-boundary diagnostics across outer folds. |

### Which observations are valid matching candidates?

Cross-validation produces one prediction for every observation, but
those predictions do not all share the same training set. A prediction
is therefore matched only against observed patterns from its own outer
test fold. Comparing it with an observation from another fold would
usually put a training target into the candidate set and can bias
identification rank below chance even for an intercept-only model.

Accordingly, `pattern_discrimination` and `pattern_rank_percentile` use
only candidates withheld together. `rdm_correlation` uses only
observation pairs withheld together. The diagonal `pattern_correlation`,
global reconstruction metrics, and per-voxel temporal correlations
retain their ordinary meanings. Permutation tests shuffle prediction
labels within outer test folds. If RDM vectors are retained, they keep
their full lower-triangle shape but mark cross-fold pairs as `NA`;
[`feature_rsa_rdm_vectors()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_rdm_vectors.md)
also returns the fold assignment.

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
#> 1      1               0.285                  0.282                   0.756
#> 2      2               0.325                  0.336                   0.789
#> 3      3               0.251                  0.282                   0.738
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
supports four estimators. The right choice depends on the feature space:

| Method | Best when | Notes |
|----|----|----|
| `"pls"` (default) | Features are correlated; you want supervised dimensionality reduction. | Uses the configured `pls` PLS numerical kernel. Fits up to `max_comps` components. |
| `"pca"` | You want to reduce features unsupervised first, then regress on PCs. | Uses the configured `pls` PCR numerical kernel (SVD-PCR by default). |
| `"ridge"` | Effects are plausibly dense; you want fast, stable L2 shrinkage for many voxel responses. | Economy-SVD multi-response ridge. Tune with `lambda_selection`; coefficients are not sparse. |
| `"glmnet"` | You expect sparse or grouped feature effects and need an elastic-net penalty. | Multivariate Gaussian elastic net via `glmnet`. Tune with `alpha`, `lambda`, `cv_glmnet`. |

## How should you select the number of components?

PLS and PCR fit up to `max_comps`; `ncomp_selection` decides how many
are used for each outer-fold prediction. This tuning happens **inside**
the outer cross-validation loop, so held-out outer observations never
choose their own component count.

| Selection | Inner work per outer training fold | What it estimates | When to use it |
|----|---:|----|----|
| `"blocked"` | One fit per training block, plus the final fit | Leave-one-block-out MSE or pattern-relative performance; each block contributes one score | Preferred when runs, sessions, subjects, or another genuine independence unit are available |
| `"loo"` (default) | One fit per training observation, plus the final fit | Leave-one-observation-out MSE | Backward-compatible choice for exchangeable, independent rows; expensive for long designs |
| `"pve"` | One fit | Cumulative predictor variance explained, stopping at `pve_threshold` (default 0.9) | A fast heuristic when held-out tuning is not required |
| `"max"` | One fit | No tuning; always use the available `max_comps` | A fixed-complexity analysis chosen in advance |

`ncomp_objective = "mse"` remains the default in this release,
preserving the existing selector and its one-standard-error rule: choose
the fewest components whose mean segment MSE is within one standard
error of the minimum. When the scientific target is relative pattern
prediction, blocked selection can instead use
`ncomp_objective = "pattern_discrimination"`. This maximizes the mean
within-block correlation advantage of the correct neural pattern over
incorrect candidate patterns. It is continuous, centered at zero under
exchangeable null matching, and scales linearly in the number of
candidates (at fixed voxel count) without forming a full candidate
correlation matrix. `"pattern_rank_percentile"` remains available when
identification rank itself is the planned estimand, but it is more
discrete and tie-sensitive.

Both pattern-relative objectives require `ncomp_selection = "blocked"`
and at least two observations in every inner validation block. They
default to the empirical optimum (`ncomp_one_se = FALSE`). Inner
predictions are converted back to the original neural-response scale
before scoring, matching the reported outer-fold metrics. Feature and
neural centering and scaling are re-estimated from each inner training
split, so the held-out block does not set preprocessing parameters. With
unequal block sizes, blocks receive equal weight rather than
observations receiving equal weight.

The default remains `"loo"` for compatibility. On blocked neuroimaging
data, pass the same observation-aligned `block_var` to
[`feature_rsa_design()`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md)
and use `ncomp_selection = "blocked"`, as in the minimal example. Every
outer training fold must retain at least two non-empty blocks; rMVPA
fails the fold rather than silently reverting to a different selector.

``` r

pca_discrimination_mod <- feature_rsa_model(
  dataset = dset,
  design = design,
  method = "pca",
  ncomp_selection = "blocked",
  ncomp_objective = "pattern_discrimination",
  ncomp_one_se = FALSE,
  crossval = cv
)
```

The implementation streams validation scores and retains only the
compact coefficient path needed for prediction. This reduces transient
memory but does not change the number of fits implied by the selector.
In particular, exact LOO remains intrinsically proportional to the
number of training observations. Parallel workers keep one target matrix
plus fold indices instead of serializing overlapping target slices for
every fold.

## When should you use fast ridge?

Use `method = "ridge"` when the scientific model is a dense linear
mapping from features to voxel responses and L2 shrinkage is acceptable.
Ridge is a distinct estimator: it is neither approximate PLS nor sparse
elastic net. It solves, after outer-fold standardisation,

``` math
n^{-1}\lVert Y-XB\rVert_F^2 + \lambda\lVert B\rVert_F^2,
```

with an unpenalised intercept. The normal equations therefore add
`n * lambda` to the feature cross-product, so a supplied lambda has the
same normalised meaning in folds of different sizes. `lambda = 0`
requests the minimum-norm unpenalised solution.

| `lambda_selection` | Inner work per outer training fold | Assumption and use |
|----|---:|----|
| `"gcv"` (default) | One economy SVD, scoring the full lambda grid | Fastest choice when rows can be treated as exchangeable. Minimises generalized cross-validation error. |
| `"loo"` | One economy SVD plus analytic PRESS scoring | Exact LOO for the ridge estimator on the fixed outer-fold scale. MSE tuning uses the one-standard-error rule by default. |
| `"blocked"` | One economy SVD per retained training block | Preferred for dependent runs, sessions, or subjects. Re-estimates centering and scaling inside every inner split. Can tune MSE, pattern discrimination, or identification rank. |
| `"fixed"` | One final fit | Use one prespecified non-negative lambda; no data-driven inner tuning. |

GCV and analytic LOO avoid a refit per observation. They still treat
individual observations as the validation unit, and LOO keeps the
outer-fold column scales fixed. If within-run dependence or fully nested
preprocessing matters, use `"blocked"`; speed is not a reason to change
the independence unit. The default grid is 49 fixed log-spaced values
from `1e-6` to `1e2`.

`lambda_objective = "mse"` remains the default in this release and is
the only objective available for GCV and analytic LOO. For relative
pattern prediction, use blocked tuning with
`lambda_objective = "pattern_discrimination"`. The score is the same
correct-minus-incorrect correlation advantage reported as
`pattern_discrimination`, is evaluated in the original neural-response
space, and uses a linear-time scorer. If identification rank itself is
the primary estimand, `"pattern_rank_percentile"` is also available.
Every candidate set is confined to one held-out inner block.
Pattern-relative tuning uses the empirical optimum by default; set
`lambda_one_se = TRUE` only if the stronger-regularisation tie-break is
part of the planned analysis.

``` r

ridge_discrimination_mod <- feature_rsa_model(
  dataset = dset,
  design = design,
  method = "ridge",
  lambda_selection = "blocked",
  lambda_objective = "pattern_discrimination",
  lambda_one_se = FALSE,
  crossval = cv
)
```

``` r

ridge_rank_mod <- feature_rsa_model(
  dataset = dset,
  design = design,
  method = "ridge",
  lambda_selection = "blocked",
  lambda_objective = "pattern_rank_percentile",
  lambda_one_se = FALSE,
  crossval = cv
)
```

``` r

ridge_mod <- feature_rsa_model(
  dataset = dset,
  design = design,
  method = "ridge",
  lambda_selection = "blocked",
  crossval = cv
)
ridge_res <- run_regional(ridge_mod, region_mask)
ridge_res$performance_table[
  , c("roinum", "pattern_correlation", "rdm_correlation",
      "median_lambda", "mean_effective_df", "mean_nonintercept_df",
      "lambda_min_boundary_fraction", "lambda_max_boundary_fraction")
]
#> # A tibble: 3 × 8
#>   roinum pattern_correlation rdm_correlation median_lambda mean_effective_df
#>    <int>               <dbl>           <dbl>         <dbl>             <dbl>
#> 1      1               0.293           0.636          3.40              1.95
#> 2      2               0.326           0.688          3.16              2.12
#> 3      3               0.211           0.554          6.81              1.65
#> # ℹ 3 more variables: mean_nonintercept_df <dbl>,
#> #   lambda_min_boundary_fraction <dbl>, lambda_max_boundary_fraction <dbl>
```

`median_lambda` is the median selected penalty across outer folds.
`mean_effective_df` is the corresponding mean effective model dimension,
including the intercept; it is the ridge analogue of a complexity
diagnostic, not a component count. `mean_nonintercept_df` subtracts that
intercept. `lambda_min_boundary_fraction` and
`lambda_max_boundary_fraction` report the fraction of outer folds
selecting an endpoint of the supplied grid. A high boundary fraction is
a diagnostic to inspect the grid and the tuning estimand, not evidence
by itself that the ridge calculation failed.

``` r

methods <- c("pls", "pca", "ridge")
summary_rows <- lapply(methods, function(m) {
  spec <- feature_rsa_model(
    dataset = dset, design = design, method = m, crossval = cv,
    ncomp_selection = "blocked", lambda_selection = "blocked"
  )
  out <- run_regional(spec, region_mask)$performance_table
  data.frame(
    method = m,
    pattern_correlation = mean(out$pattern_correlation),
    rdm_correlation = mean(out$rdm_correlation),
    mse = mean(out$mse)
  )
})
do.call(rbind, summary_rows)
#>   method pattern_correlation rdm_correlation      mse
#> 1    pls           0.2870185       0.5951209 1.229655
#> 2    pca           0.2855695       0.6046905 1.232471
#> 3  ridge           0.2765517       0.6258379 1.250525
```

## What did the ridge benchmark show?

The checked-in characterization uses 20 independently seeded dense
linear fixtures. Each has 120 observations, 16 correlated features, 100
voxel responses, six blocks, and a never-tuned-on final held-out block.
All native thread limits were one. The prospective accuracy margin
required median MSE no more than 5% worse than PLS LOO, 90th-percentile
MSE no more than 20% worse, and pattern/RDM correlation losses no larger
than .01 at the median or .03 at the lower decile.

| Estimator / selector | Median ms | Speedup vs PLS LOO | MSE ratio | Pattern delta | RDM delta |
|:---|---:|---:|---:|---:|---:|
| PLS LOO | 90.00 | 1.0 | 1.000 | 0.0000 | 0.0000 |
| PLS blocked | 8.00 | 11.1 | 1.000 | 0.0000 | 0.0000 |
| Ridge GCV | 2.17 | 41.3 | 0.877 | 0.0022 | 0.0026 |
| Ridge analytic LOO | 8.17 | 11.1 | 0.887 | 0.0018 | 0.0021 |
| Ridge blocked | 8.67 | 10.3 | 0.888 | 0.0018 | 0.0018 |
| glmnet ridge blocked CV | 676.00 | 0.1 | 2.624 | -0.0100 | -0.0126 |

Twenty-seed dense-linear characterization on Apple M3 Max. Accuracy
columns are relative to PLS LOO; timings are descriptive. {.table}

All three ridge selectors cleared the declared margin. Median MSE was
11–12% lower than PLS LOO, while median pattern and RDM correlations
were slightly higher. GCV was about 39 times faster than PLS LOO;
analytic LOO and blocked ridge were about 10 and 9 times faster. Blocked
PLS was also about 11 times faster than PLS LOO and had essentially
identical accuracy here, so users who specifically need PLS should
prefer a scientifically valid block unit rather than changing estimators
only for speed.

The `glmnet` result is not a general ranking: it reports that package’s
blocked-CV choice on this particular dense fixture, where it selected a
different penalty and missed the fixture-specific equivalence margin.
Ridge does not replace elastic net when sparsity is part of the
estimand. Likewise, these simulations do not prove equivalence on a new
dataset; compare outer-fold metrics on your own ROIs before changing a
preregistered method.

Run `inst/benchmarks/feature_rsa_hotpaths.R` from a source checkout to
regenerate the receipt. The driver validates numerical parity before
timing, records a source fingerprint and dirty-worktree status, and
includes the former materialized blocked-ridge path as an oracle. The
implementation uses base R’s optimized SVD and matrix multiplication;
profiling did not justify an Rcpp layer around those same kernels.

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
