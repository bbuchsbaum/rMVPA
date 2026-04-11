# Feature-sets design (mvpa_design extension for continuous regression)

\`feature_sets_design()\` defines a small \`mvpa_design\` extension that
attaches grouped continuous predictors for stimulus-\>brain regression.

## Usage

``` r
feature_sets_design(
  X_train,
  X_test = NULL,
  block_var_test = NULL,
  target_builder = NULL,
  target_builder_data = NULL,
  n_test = NULL,
  split_by = NULL
)
```

## Arguments

- X_train:

  A \`feature_sets\` object (train/source predictors).

- X_test:

  Optional \`feature_sets\` object (test/target predictors).

- block_var_test:

  Optional integer/factor vector of length T_rec defining test/target
  blocks (typically runs). Used by models that do test-time CV.

- target_builder:

  Optional function that rebuilds target-domain predictors per outer
  target fold. It is called with named arguments including \`X_train\`,
  \`X_test\`, \`train_idx\`, \`test_idx\`, \`fold_id\`,
  \`block_var_test\`, \`n_test\`, \`builder_data\`, and \`design\`;
  callbacks may accept only the subset they need. The return value must
  resolve to target predictors with the same feature/set layout as
  \`X_train\`.

- target_builder_data:

  Optional object passed through to \`target_builder\` as
  \`builder_data\`.

- n_test:

  Optional target-domain row count used when \`X_test\` is not supplied.
  If omitted, it is inferred from \`X_test\` or \`block_var_test\`.

- split_by:

  Optional split definition passed to \`mvpa_design\`.

## Value

An object inheriting from \`mvpa_design\` with class
\`feature_sets_design\`.

## Details

This mirrors the pattern used by the archived \`hrfdecoder\` design API:
rMVPA's core infrastructure expects \`mvpa_design\` fields like
\`train_design\`, \`y_train\`, etc. For continuous regression with
external test domains, those fields are largely bookkeeping, while the
*actual* predictors are carried explicitly as \`feature_sets\` objects.

**Where the data live.**

- Predictors live on the design: \`design\$X_train\` and
  \`design\$X_test\` (both \`feature_sets\`).

- Responses live on the dataset: \`dataset\$train_data\` and
  \`dataset\$test_data\` (TRxvoxel matrices per ROI).

**cv_labels semantics.** The returned object passes \`cv_labels =
1:T_enc\` to \`mvpa_design()\`. These integer indices are used only for
length validation and fold-construction bookkeeping in the generic rMVPA
iterators; they are not meaningful training targets. The actual
predictors are carried as \`design\$X_train\` / \`design\$X_test\`, and
models built for this design (e.g. \`grouped_ridge_da_model()\` /
\`banded_ridge_da_model()\`) retrieve them from those fields rather than
from \`cv_labels\`.

**Recall blocks.** \`block_var_test\` is stored for convenience as
\`design\$block_var_test\`. Models can use it to define blocked
cross-validation on test/target time (e.g. leave-one-run-out), and fall
back to contiguous folds when only a single run is present.

**Fold-aware target builders.** For unbiased domain-adaptation
workflows, \`feature_sets_design()\` can store a \`target_builder\`
callback instead of a single fixed \`X_test\`. The callback is invoked
separately for each outer target fold and receives the source predictors
plus the target train/test row indices. It must return target-side
predictors in the original target row order, either as a
\`feature_sets\` object, a numeric matrix, or a list containing
\`gamma\`, \`X\`, or \`X_test\`. This allows upstream matching/alignment
to be re-fit on target-train rows only, while held-out target rows
remain untouched until evaluation.

## See also

[`feature_sets`](http://bbuchsbaum.github.io/rMVPA/reference/feature_sets.md),
[`expected_features`](http://bbuchsbaum.github.io/rMVPA/reference/expected_features.md),
[`banded_ridge_da_model`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md),
[`grouped_ridge_da_model`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da_model.md),
[`banded_ridge_da`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md),
[`grouped_ridge_da`](http://bbuchsbaum.github.io/rMVPA/reference/banded_ridge_da.md)

## Examples

``` r
# Train predictors (TR x features), split into named sets:
X_enc <- matrix(rnorm(20 * 8), 20, 8)
fs_enc <- feature_sets(X_enc, blocks(low = 3, sem = 5))

# Test predictors (TR x features), for example from a soft alignment:
gamma <- matrix(runif(10 * 20), 10, 20)
gamma <- gamma / rowSums(gamma)
fs_rec <- expected_features(fs_enc, gamma, drop_null = FALSE, renormalize = TRUE)

des <- feature_sets_design(fs_enc, fs_rec, block_var_test = rep(1:2, each = 5))
```
