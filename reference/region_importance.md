# Region Importance via Random Subset Comparison

Assess each feature's (region/parcel/voxel) contribution to model
performance by comparing cross-validated accuracy when the feature is
included vs. excluded across many random feature subsets.

## Usage

``` r
region_importance(model_spec, ...)

# S3 method for class 'mvpa_model'
region_importance(
  model_spec,
  n_iter = 200,
  subset_fraction = 0.5,
  metric = NULL,
  aggregation = c("mean", "sum", "maxabs"),
  ...
)
```

## Arguments

- model_spec:

  An `mvpa_model` specification.

- ...:

  Additional arguments (currently unused).

- n_iter:

  Number of random subset iterations (default 200).

- subset_fraction:

  Fraction of features sampled per iteration (default 0.5).

- metric:

  Character name of performance metric to use (e.g. "Accuracy"). NULL
  uses the first metric returned by `compute_performance`.

- aggregation:

  How to aggregate multi-basis feature importance (default "mean").

## Value

A `region_importance_result` object.

A `region_importance_result` object.

## Details

**Interpretation.** `region_importance` measures each feature's
*marginal predictive contribution*: how much cross-validated performance
improves, on average, when the feature is included versus excluded
across random feature subsets. This is an approximation of Shapley
values and is purely a **decoding** (backward) measure – it never
examines model weights.

**Haufe et al. (2014) considerations.** Because this method works
through out-of-sample performance deltas rather than weight
interpretation, it sidesteps the core failure mode described by Haufe et
al. (2014), where backward-model weights are misinterpreted as
activation patterns.

However, *suppressor variables* – features that improve classification
by cancelling correlated noise rather than carrying signal – will
receive positive importance, since they genuinely improve
generalization. This is correct from a decoding perspective ("which
features help classify?") but potentially misleading from a neuroscience
perspective ("where does the signal originate?").

**Complementary methods.** For forward-model (activation-pattern)
interpretation of linear models, use
[`haufe_importance`](http://bbuchsbaum.github.io/rMVPA/reference/haufe_importance.md)
directly or rely on the `importance_vector` returned by
[`run_global`](http://bbuchsbaum.github.io/rMVPA/reference/run_global.md),
which uses
[`model_importance`](http://bbuchsbaum.github.io/rMVPA/reference/model_importance.md)
internally. For non-linear models such as random forests,
`region_importance` is the recommended approach.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(5,5,5), 40, nlevels=2, blocks=3)
  cval <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("sda_notune")
  mspec <- mvpa_model(mdl, ds$dataset, ds$design,
    "classification", crossval=cval)
  imp <- region_importance(mspec, n_regions=5, n_subsets=10)
# }
```
