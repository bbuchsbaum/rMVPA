# ERA-RSA: Encoding-Retrieval Similarity and ER Geometry

Combines first-order encoding-retrieval similarity (ERA) with
second-order RSA between encoding and retrieval representational
geometries. Works with
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
and
[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
using the standard rMVPA iterators.

## Usage

``` r
era_rsa_model(
  dataset,
  design,
  key_var,
  phase_var,
  encoding_level = NULL,
  retrieval_level = NULL,
  distfun = cordist(method = "pearson"),
  rsa_simfun = c("pearson", "spearman"),
  confound_rdms = NULL,
  partial_against = "run",
  include_diag = TRUE,
  item_block = NULL,
  item_lag = NULL,
  item_run_enc = NULL,
  item_run_ret = NULL,
  global_nuisance = FALSE,
  require_run_metadata = FALSE,
  ...
)
```

## Arguments

- dataset:

  An mvpa_dataset with train_data (encoding) and test_data (retrieval).

- design:

  An mvpa_design describing trial structure with train/test designs.

- key_var:

  Column name or formula giving the item key that links encoding and
  retrieval trials (e.g., ~ ImageID).

- phase_var:

  Column name or formula giving phase labels (must include encoding and
  retrieval levels if using a single-phase dataset; for the default
  two-dataset usage, this is still parsed for consistency but not
  required for operations).

- encoding_level:

  Level of phase_var to treat as encoding (default: first level).

- retrieval_level:

  Level of phase_var to treat as retrieval (default: second level).

- distfun:

  A distfun used to compute within-phase RDMs (e.g.,
  cordist("pearson")).

- rsa_simfun:

  Character: similarity for comparing RDMs, one of "pearson" or
  "spearman".

- confound_rdms:

  Optional named list of KxK matrices or "dist" objects encoding
  item-level nuisance/model RDMs (e.g., block, run, time). Rows and
  columns should correspond to item keys (levels of `key_var`). When
  `run_enc` and `run_ret` entries are present they are used to compute
  `geom_cor_run_partial`, the ER geometry correlation after regressing
  out these run confounds. The same list also supplies candidate
  nuisance RDMs for `geom_cor_partial`, selected by `partial_against`.

- partial_against:

  Character vector selecting nuisance groups or exact `confound_rdms`
  names used for the general `geom_cor_partial` metric. Recognized
  groups are `"run"`, `"time"`, `"block"`, `"category"`, `"global"`, and
  `"all"`. Exact confound names such as `"time_enc"` are also allowed.
  The default `"run"` preserves the legacy run-partial interpretation
  while exposing the result under the more general `geom_cor_partial`
  name. If `global_nuisance` is enabled and `partial_against` is not
  supplied explicitly, the effective default is `c("run", "global")`.

- include_diag:

  Logical retained for API compatibility. Diagonal ERA metrics are
  always retained, and the off-diagonal mean used by
  `era_diag_minus_off` always excludes matching-item diagonal entries.

- item_block:

  Optional factor of per-item blocks, aligned to item keys. Typically
  derived by aggregating a trial-level block/run variable in
  `design$train_design` to the item level (e.g., modal block per item).
  Used to compute block-limited ERA contrasts
  (`era_diag_minus_off_same_block` / `era_diag_minus_off_diff_block`).

- item_lag:

  Optional numeric per-item retrieval-minus-encoding lag, aligned to
  item keys. Often derived from trial onsets (e.g., mean retrieval onset
  minus mean encoding onset per item). Used to compute `era_lag_cor`,
  the correlation between ERA diagonal and lag.

- item_run_enc:

  Optional factor of per-item encoding runs, aligned to item keys (e.g.,
  modal run for encoding trials of each item). Combined with
  `item_run_ret` to compute `geom_cor_xrun`, ER geometry restricted to
  item pairs differing in both encoding and retrieval run.

- item_run_ret:

  Optional factor of per-item retrieval runs, aligned to item keys. See
  `item_run_enc` for how it is used.

- global_nuisance:

  Logical or pre-supplied list controlling whole-mask global similarity
  nuisance. `FALSE` (default) disables it. `TRUE` computes K x K
  item-level RDMs (`global_enc`, `global_ret`) over the full
  `dataset$mask` once at construction time and adds them to
  `confound_rdms`. They are then picked up by `geom_cor_partial` when
  `partial_against` includes `"global"` (or `"all"`). A pre-computed
  list with elements `D_enc`/`enc` and `D_ret`/`ret` can be supplied
  directly for non-standard dataset backends. Caveat: each ROI/sphere is
  part of the global mask, so for large regional ROIs covering most of
  the mask the residualization partially removes the local signal.

- require_run_metadata:

  Logical; if `TRUE`, missing item-level run or block metadata becomes
  an error rather than a warning. Use this when a downstream analysis
  depends on `era_diag_minus_off_same_block`,
  `era_diag_minus_off_diff_block`, or `geom_cor_xrun` — the constructor
  will refuse to silently produce schemas where those metrics are
  guaranteed to be `NA`. Default `FALSE` (warn only).

- ...:

  Additional fields stored on the model spec.

## Value

A model spec of class `"era_rsa_model"` compatible with
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
and
[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md).
When fit, the spec emits the scalar metrics documented in **Metrics**.

## Details

Key outputs per ROI/searchlight sphere include: - First-order ERA: top-1
accuracy, diagonal mean, diagonal-minus-off. - Second-order geometry:
correlation between encoding and retrieval RDMs. - Optional
confound-aware metrics and diagnostics (block/lag/run).

`era_rsa_model()` itself returns a model specification. The scalar
outputs below are emitted when that specification is evaluated with
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
or
[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md):
regional analyses place them in `result$performance_table`; searchlight
analyses expose them as metric maps listed in `result$metrics`.

## Metrics

For each ROI / searchlight center, `era_rsa_model` emits a set of scalar
metrics that are turned into spatial maps by
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
and
[`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md):

- n_items:

  Number of unique item keys contributing to this ROI/sphere (i.e.,
  length of the common encoding-retrieval item set).

- era_top1_acc:

  Top-1 encoding\\\rightarrow\\retrieval accuracy at the item level:
  fraction of retrieval trials whose most similar encoding pattern (over
  items) has the same `key_var`.

- era_diag_mean:

  Mean encoding-retrieval similarity for matching items (mean of the
  diagonal of the encoding\\\times\\retrieval similarity matrix).

- era_diag_minus_off:

  Diagonal-minus-off-diagonal ERA contrast: `era_diag_mean` minus the
  mean similarity to all non-matching items, capturing how much
  same-item pairs stand out from other pairs.

- geom_cor:

  Correlation between encoding and retrieval representational
  geometries: correlation (Pearson or Spearman, per `rsa_simfun`)
  between the vectorised lower triangles of the encoding and retrieval
  RDMs.

- geom_cor_partial:

  General nuisance-partial ER geometry correlation. This residualizes
  both encoding and retrieval geometry vectors against the confounds
  selected by `partial_against`, then correlates the residuals. Run
  confounds can be supplied either as `confound_rdms$run_enc`/`run_ret`
  or derived from `item_run_enc`/`item_run_ret`.

- era_diag_minus_off_same_block:

  Block-limited ERA contrast when `item_block` is supplied: diagonal ERA
  minus the mean similarity to other items in the same block (e.g.,
  run/condition), averaged over items.

- era_diag_minus_off_diff_block:

  Cross-block ERA contrast when `item_block` is supplied: diagonal ERA
  minus the mean similarity to items in different blocks.

- era_lag_cor:

  Lag-ERA correlation when `item_lag` is supplied: Spearman correlation
  between diagonal ERA values and the per-item lag (e.g., retrieval
  minus encoding onset), using complete cases only.

- geom_cor_run_partial:

  Run-partial ER geometry correlation, when run-level confounds are
  supplied via `confound_rdms$run_enc`/`run_ret` or derivable from
  `item_run_enc`/`item_run_ret`. Computed as the correlation between
  encoding and retrieval RDMs after regressing out those run RDMs.

- geom_cor_xrun:

  Cross-run-only ER geometry correlation, when `item_run_enc` and
  `item_run_ret` are supplied: correlation between encoding and
  retrieval RDMs restricted to item pairs that differ in both encoding
  and retrieval run.

- beta\_\*:

  When `confound_rdms` are provided, additional `beta_<name>` terms give
  the regression coefficients from a linear model predicting retrieval
  geometry from encoding geometry and the confound RDMs (one coefficient
  per nuisance/model RDM).

- sp\_\*:

  If `run_lm_semipartial()` is available, `sp_<name>` terms provide
  semi-partial R\\^2\\-like diagnostics for each confound RDM,
  quantifying unique variance explained in retrieval geometry.

## Trial-level vs. item-level metadata

`block_var` on
[`mvpa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
is *trial-level* metadata (one entry per row of the design table) and is
used by cross-validation, not by the ERA item-level metrics. The
block/run-aware metrics (`era_diag_minus_off_same_block`,
`era_diag_minus_off_diff_block`, `geom_cor_xrun`) are computed from
*item-level* vectors indexed by levels of `key_var`: `item_block`,
`item_run_enc`, and `item_run_ret`. These must be supplied directly
here; passing `block_var = ~run` to
[`mvpa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
alone will not enable them, and the model will warn that those metrics
will be `NA`. See
[`era_rsa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/era_rsa_design.md)
for a helper that builds these vectors from the design table.

For external train/test phases (encoding vs. retrieval), run labels
often collide across phases (both phases may have runs `1, 2, 3` that
correspond to different scans). When `item_run_enc` and `item_run_ret`
share atomic values that are not phase-scoped, supply phase-prefixed
labels such as `enc_1` / `ret_1` so cross-phase equality tests are not
spuriously satisfied.

## Examples

``` r
if (FALSE) { # \dontrun{
  # See vignette for complete ERA-RSA workflow
} # }
```
