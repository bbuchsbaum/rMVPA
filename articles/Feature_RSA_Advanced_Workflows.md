# Feature-RSA Advanced Workflows

Once ordinary feature-RSA is working, the next problem is usually
transfer. Sometimes you want to stay within one ROI and ask whether a
feature space explains held-out data. Sometimes you want to move across
cognitive states. Sometimes you want to ask whether geometry learned in
one ROI appears in another.

Those are different estimands. This vignette is a short guide to
choosing the right one and interpreting the result without mixing them
together.

## What is the shortest map of the workflow space?

``` r
advanced_summary <- build_advanced_summary()

knitr::kable(advanced_summary, digits = 3)
```

| workflow                 | core_function                                            | summary_value |
|:-------------------------|:---------------------------------------------------------|--------------:|
| within_roi               | feature_rsa_model()                                      |         0.036 |
| cross_state_da           | feature_rsa_da_model()                                   |         0.090 |
| cross_roi_generalization | feature_rsa_cross_connectivity()                         |         0.000 |
| offset_diagnostics       | feature_rsa_cross_connectivity(adjust = ‘double_center’) |         0.162 |

The table gives you four practical signals:

- a within-ROI held-out geometry score
- the gain from adapting across states instead of using source rows
  alone
- the diagonal advantage in the ROI_i -\> ROI_j matrix
- the spread of source and target offsets that can create stripe
  patterns

## When is standard feature-RSA enough?

Stay with
[`feature_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
when your question is local to one ROI and one state. You fit the model
in the usual way, hold out data within that state, and read the ROI-wise
geometry scores from
[`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md).

``` r
base_fit <- build_base_feature_rsa()

base_fit$result$performance_table[, c("roinum", "rdm_correlation")]
#> # A tibble: 4 × 2
#>   roinum rdm_correlation
#>    <int>           <dbl>
#> 1      1        -0.00825
#> 2      2         0.0563 
#> 3      3         0.00301
#> 4      4         0.0912
```

This is the right default when you want to know whether a feature space
explains one parcel on its own terms.

## When do you need domain adaptation?

Switch to
[`feature_rsa_da_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_da_model.md)
when the feature space is shared but the state changes. The canonical
example is encoding to recall: the rows still refer to the same latent
content, but the source and target mappings need not be identical.

``` r
da_perf <- build_small_da_example()

knitr::kable(da_perf, digits = 3)
```

| model       | target_rdm_correlation | target_r2_full |
|:------------|-----------------------:|---------------:|
| source_only |                  0.857 |          0.566 |
| coupled_da  |                  0.948 |          0.835 |

The important comparison is not just target prediction in the abstract.
It is whether a model that can adapt on target-train rows recovers more
held-out target geometry than a source-only baseline. When the
target-side representation is itself estimated from target data, use a
fold-aware `target_builder` as shown in
[`vignette("Feature_RSA_Domain_Adaptation")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Domain_Adaptation.md).

## When do you want ROI_i -\> ROI_j transfer?

Use
[`feature_rsa_cross_connectivity()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_cross_connectivity.md)
when you want a full asymmetric ROI x ROI matrix. The rows are source
ROIs, the columns are target ROIs, and the off-diagonals tell you
whether geometry learned in one parcel matches held-out observed
geometry in another.

``` r
knitr::kable(round(base_fit$cross_raw, 2))
```

|     1 |    2 |     3 |    4 |
|------:|-----:|------:|-----:|
| -0.01 | 0.06 | -0.01 | 0.13 |
|  0.01 | 0.06 |  0.03 | 0.03 |
| -0.03 | 0.09 |  0.00 | 0.07 |
|  0.01 | 0.03 |  0.00 | 0.09 |

This is the right analysis when ROI size differs and voxel-space weight
transfer would be ill-posed, but you still want to ask whether
representational geometry generalizes across parcels.

## What do you do about source-dominant and target-dominant ROIs?

If the ROI x ROI matrix shows vertical or horizontal striping, some ROIs
are globally easy targets or globally strong sources. That can be a real
feature of the data, but it can also swamp the pair-specific structure
you actually care about.

``` r
offset_table <- build_offset_example()

knitr::kable(offset_table, digits = 2)
```

| roinum | source_offset | target_offset |
|-------:|--------------:|--------------:|
|      1 |          0.13 |          0.13 |
|      2 |          0.08 |          0.11 |
|      3 |          0.00 |          0.05 |
|      4 |         -0.21 |         -0.29 |

Use `adjust = "double_center"` when you want to remove additive source
and target main effects from the ROI x ROI matrix itself. Use
`adjust = "residualize_mean"` when you think the stripes are driven by
one dominant common geometry that should be removed before the
cross-correlation is computed.

## Decision rules

- Use
  [`feature_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md)
  for one ROI, one state, held-out evaluation in the same domain.
- Use
  [`feature_rsa_da_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_da_model.md)
  when the feature space is shared across states and the target mapping
  may shift.
- Use
  [`feature_rsa_cross_connectivity()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_cross_connectivity.md)
  when you want ROI-to-ROI generalization at the level of predicted
  versus observed geometry.
- Use
  [`feature_rsa_connectivity()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_connectivity.md)
  when you want a symmetric similarity matrix of predicted geometries
  themselves.
- Use `adjust = "double_center"` or `adjust = "residualize_mean"` only
  after you are clear about which nuisance pattern you are trying to
  remove.

## Next steps

For the full cross-state workflow, including fold-aware target
rebuilding, see
[`vignette("Feature_RSA_Domain_Adaptation")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Domain_Adaptation.md).
For the full ROI-to-ROI workflow and stripe-mitigation examples, see
[`vignette("Feature_RSA_Connectivity")`](http://bbuchsbaum.github.io/rMVPA/articles/Feature_RSA_Connectivity.md).
