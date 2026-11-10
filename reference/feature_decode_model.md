# Feature decoding model: predict a feature space from activity patterns

The pattern-level complement of \[banded_ridge_model()\]. Within each
region or searchlight, the feature columns of a
\[feature_sets_design()\] are predicted from the local activity pattern
by cross-validated ridge regression, and the held-out predictions are
summarised over the declared structure of the feature space: every set,
every family of sets, and every contrast (see \[feature_sets()\]).
Family and contrast summaries are computed from the same held-out
predictions, so they need no additional fits.

## Usage

``` r
feature_decode_model(
  dataset,
  design,
  outer_crossval = NULL,
  tune_crossval = NULL,
  alphas = 10^seq(-1, 5),
  purge = 0L,
  nuisance = NULL,
  nuisance_fit = c("train", "within"),
  partial = FALSE,
  n_perm = 0L,
  perm_strategy = c("circular_shift", "block_shuffle"),
  perm_block = NULL,
  min_shift = 1L,
  min_voxels = 2L,
  seed = 1L,
  ...
)
```

## Arguments

- dataset:

  An \`mvpa_dataset\`; rows of its training data are the observations
  (for example volumes of one run).

- design:

  A \[feature_sets_design()\] whose training features are the decoding
  targets, with the same number of rows as the dataset.

- outer_crossval, tune_crossval, purge:

  Cross-validation, as in \[banded_ridge_model()\]: an integer fold
  count, a \`cross_validation\` object, or explicit folds;
  \`tune_crossval\` sets the inner folds used to choose the ridge
  penalty; \`purge\` drops training rows within that many rows of a
  held-out row.

- alphas:

  Candidate ridge penalties. Each feature column selects its own penalty
  inside each outer-training set.

- nuisance, nuisance_fit:

  Optional nuisance matrix (rows = observations) regressed out of
  patterns and targets inside every split, and how held-out rows are
  residualized; see \[banded_ridge_model()\].

- partial:

  If \`TRUE\`, also report partial accuracies (see Metrics).

- n_perm:

  Number of time-shift permutations; \`0\` (default) disables the null.

- perm_strategy, perm_block, min_shift:

  Permutation scheme: \`"circular_shift"\` (default) or
  \`"block_shuffle"\`, the block length for the latter, and the smallest
  circular shift in rows.

- min_voxels:

  Regions with fewer voxels are returned as errors.

- seed:

  Integer seed for fold construction and permutations.

- ...:

  Passed to the model specification.

## Value

A model specification of class \`feature_decode_model\`.

## Details

The model shares its cross-validation, temporal purge, nuisance handling
and ridge solvers with \[banded_ridge_model()\], with the roles of the
two matrices exchanged. It implements \[fit_roi()\], so
\[run_searchlight()\] and \[run_regional()\] work unchanged.

## Metrics

One value per region or searchlight centre:

- \`all_r2\`, \`all_cor\`:

  Mean held-out R2 and correlation over all feature columns.

- \`\<set\>\_r2\`, \`\<set\>\_cor\`:

  The same, over the columns of one set or family. Every feature column
  counts equally.

- \`\<contrast\>\_r2\`, \`\<contrast\>\_cor\`:

  Accuracy for the contrast series, i.e. observed and predicted features
  projected onto the contrast.

- \`\<set\>\_partial_r2\`, \`\<set\>\_partial_cor\`:

  With \`partial = TRUE\`: accuracy for a set or family after every
  other feature column has been regressed out of both the pattern and
  the targets, with coefficients fitted on each training split. This is
  the pattern-level analogue of unique variance.

- \`\<metric\>\_p\`:

  With \`n_perm \> 0\`: one-sided p-value against a null in which the
  predictions are shifted in time against the observations.

## See also

\[banded_ridge_model()\], \[feature_sets()\], \[run_searchlight()\]

## Examples

``` r
# \donttest{
set.seed(1)
n <- 80
ds <- gen_sample_dataset(c(5, 5, 5), nobs = n, blocks = 4)
X <- matrix(rnorm(n * 3), n, 3, dimnames = list(NULL, c("m1", "p1", "ac")))
fs <- feature_sets(X, by_set(c("manner", "place", "acoustic")),
                   families = list(phonology = c("manner", "place")))
des <- feature_sets_design(fs, block_var_train = rep(1:4, each = 20),
                           time_series = TRUE)
model <- feature_decode_model(ds$dataset, des, outer_crossval = 4)
res <- run_searchlight(model, radius = 4, method = "standard")
#> INFO [2026-11-10 00:19:43] Running standard searchlight with radius = 4
#> INFO [2026-11-10 00:19:43] shard backend: preparing shared memory for dataset (mvpa_image_dataset, mvpa_dataset, list)
#> INFO [2026-11-10 00:19:43] shard backend [volumetric]: shared 80 x 125 matrix (125 masked voxels)
#> INFO [2026-11-10 00:19:43] creating standard searchlight
#> INFO [2026-11-10 00:19:43] running standard searchlight iterator
#> INFO [2026-11-10 00:19:43] Using automatic searchlight batch size 125 for 125 centers (memory budget 512.0 MiB).
#> INFO [2026-11-10 00:20:02] 
#> MVPA Iteration Complete
#> - Total ROIs: 125
#> - Processed: 125
#> - Skipped: 0
#> INFO [2026-11-10 00:20:02] searchlight (standard): 125 ROIs processed (success=125, errors=0)
# }
```
