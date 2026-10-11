# Representational connectivity model (ReNA-RC)

For each ROI or searchlight, computes representational connectivity
between the ROI RDM and a seed RDM, optionally controlling for confound
RDMs (block, lag, behavior, etc.).

## Usage

``` r
repnet_model(
  dataset,
  design,
  repnet_des,
  distfun = cordist(method = "pearson"),
  simfun = c("pearson", "spearman"),
  ...
)
```

## Arguments

- dataset:

  An `mvpa_dataset`.

- design:

  An `mvpa_design`.

- repnet_des:

  Output of
  [`repnet_design`](https://bbuchsbaum.github.io/rMVPA/reference/repnet_design.md).

- distfun:

  Distance function for ROI RDM (e.g. `cordist(method = "pearson")`).

- simfun:

  Character: similarity metric between ROI and seed RDMs (`"pearson"` or
  `"spearman"`).

- ...:

  Extra fields stored on the model spec.

## Value

A model spec of class `"repnet_model"` compatible with
[`run_regional()`](https://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
and
[`run_searchlight()`](https://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md).

## Details

Per-ROI metrics: `n_items`, `n_pairs`, `conn_raw` (the ROI-seed RDM
similarity), and from the multiple regression of the ROI RDM on the seed
and confound RDMs, `beta_<term>` (regression coefficients) and
`sp_<term>` (semi-partial correlations). With confound RDMs,
`seed_conf_r2` is the R-squared of the seed RDM regressed on the
confound RDMs: `beta_seed` and `sp_seed` describe only the
`1 - seed_conf_r2` share of the seed that the confounds do not explain.

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 16, blocks = 4, nlevels = 2)
# Seed RDM over trials (e.g. from a seed region or a behavioural model)
key_ids <- ds$design$train_design$.rownum
seed_rdm <- as.matrix(dist(key_ids))
rownames(seed_rdm) <- colnames(seed_rdm) <- as.character(key_ids)
rn_des <- repnet_design(ds$design, key_var = ~ .rownum, seed_rdm = seed_rdm)
model <- repnet_model(ds$dataset, ds$design, rn_des)
```
