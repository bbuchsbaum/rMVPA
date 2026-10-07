# Representational mediation model (ReNA-RM)

For each ROI/searchlight, tests whether the ROI RDM mediates the
relationship between a predictor RDM (X) and an outcome RDM (Y),
optionally including confound RDMs.

## Usage

``` r
repmed_model(
  dataset,
  design,
  repmed_des,
  key_var,
  distfun = cordist(method = "pearson"),
  ...
)
```

## Arguments

- dataset:

  An `mvpa_dataset`.

- design:

  An `mvpa_design`.

- repmed_des:

  Output of
  [`repmed_design`](https://bbuchsbaum.github.io/rMVPA/reference/repmed_design.md).

- key_var:

  Column or formula giving item identity in the mvpa_design (e.g.
  `~ ImageID`).

- distfun:

  Distance function to compute the ROI RDM (e.g.
  `cordist(method = "pearson")` ).

- ...:

  Extra fields stored on the model spec.

## Value

A model spec of class `"repmed_model"`.

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 12, blocks = 3, nlevels = 2)
items <- as.character(sort(unique(ds$design$train_design$.rownum)))
# Predictor (X) and outcome (Y) RDMs over items
X <- as.matrix(dist(matrix(rnorm(length(items) * 2), ncol = 2)))
Y <- as.matrix(dist(matrix(rnorm(length(items) * 2), ncol = 2)))
rownames(X) <- colnames(X) <- rownames(Y) <- colnames(Y) <- items
repmed_des <- repmed_design(items = items, X_rdm = X, Y_rdm = Y)
model <- repmed_model(ds$dataset, ds$design, repmed_des, key_var = ~ .rownum)
res <- run_regional(model, ds$dataset$mask)
#> INFO [2026-10-07 05:02:23] 
#> MVPA Iteration Complete
#> - Total ROIs: 1
#> - Processed: 1
#> - Skipped: 0
#> INFO [2026-10-07 05:02:23] run_regional: 1 ROIs processed (success=1, errors=0)
res$performance_table
#> # A tibble: 1 × 10
#>   roinum n_items n_pairs   med_a med_b med_cprime  med_c med_indirect
#>    <int>   <dbl>   <dbl>   <dbl> <dbl>      <dbl>  <dbl>        <dbl>
#> 1      1      12      66 0.00669 0.231     -0.179 -0.177      0.00154
#> # ℹ 2 more variables: med_sobel_z <dbl>, med_sobel_p <dbl>
```
