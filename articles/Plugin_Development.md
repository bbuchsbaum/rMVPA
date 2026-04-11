# Building Plugin Analyses for rMVPA

``` r
library(rMVPA)
library(neuroim2)
future::plan(future::sequential)
```

## Why plugin analyses?

`rMVPA` already provides execution harnesses for searchlight and
regional analyses (ROI extraction, iteration, parallel orchestration,
and spatial map assembly). Plugin analyses let you add new method
families in a separate package while reusing that harness.

The plugin contract is small:

1.  Constructor: returns a model specification via
    [`create_model_spec()`](http://bbuchsbaum.github.io/rMVPA/reference/create_model_spec.md).
2.  `output_schema.<your_class>()`: declares metric names and layout.
3.  `fit_roi.<your_class>()`: computes one ROI result.
4.  Optional `strip_dataset.<your_class>()`: trims large fields for
    workers.

## Classifier Registry vs Analysis Plugins

These are different extension paths:

1.  Classifier registry
    ([`register_mvpa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/register_mvpa_model.md)):
    add a new estimator backend (fit/predict/prob) to use inside the
    standard `mvpa_model` analysis.
2.  Analysis plugin contract
    ([`create_model_spec()`](http://bbuchsbaum.github.io/rMVPA/reference/create_model_spec.md) +
    `fit_roi.<class>`): add a new per-ROI analysis type with custom
    logic.

Use the classifier registry when your workflow is still standard
classification/regression and only the estimator changes. Use analysis
plugins when the ROI computation itself is novel.

## Step 1: Define a constructor

``` r
toy_plugin_model <- function(dataset, design, center_stat = c("mean", "median")) {
  center_stat <- match.arg(center_stat)
  create_model_spec(
    name = "toy_plugin_model",
    dataset = dataset,
    design = design,
    center_stat = center_stat,
    compute_performance = TRUE,
    return_predictions = FALSE
  )
}
```

## Step 2: Declare the metric schema

``` r
output_schema.toy_plugin_model <- function(model) {
  list(
    score = "scalar",
    spread = "scalar"
  )
}
```

## Step 3: Implement `fit_roi`

``` r
fit_roi.toy_plugin_model <- function(model, roi_data, context, ...) {
  x <- roi_data$train_data
  score <- if (identical(model$center_stat, "median")) stats::median(x) else mean(x)
  spread <- stats::sd(as.numeric(x))

  roi_result(
    metrics = c(score = score, spread = spread),
    indices = roi_data$indices,
    id = context$id
  )
}
```

For classification/regression-style plugins, you can delegate fold
execution to
[`cv_evaluate_roi()`](http://bbuchsbaum.github.io/rMVPA/reference/cv_evaluate_roi.md)
instead of reimplementing CV loops in each plugin:

``` r
fit_roi.my_classif_plugin <- function(model, roi_data, context, ...) {
  cv_evaluate_roi(
    model_spec = model,
    roi_data = roi_data,
    context = context,
    mode = "auto"
  )
}
```

## Step 4: Register methods in-session

For a package, these methods are discovered automatically when your
package is loaded. In an interactive session, register them explicitly:

``` r
registerS3method(
  "fit_roi",
  "toy_plugin_model",
  fit_roi.toy_plugin_model,
  envir = asNamespace("rMVPA")
)
registerS3method(
  "output_schema",
  "toy_plugin_model",
  output_schema.toy_plugin_model,
  envir = asNamespace("rMVPA")
)
```

## Fast plugin validation without full pipelines

Use
[`mock_roi_data()`](http://bbuchsbaum.github.io/rMVPA/reference/mock_roi_data.md),
[`mock_context()`](http://bbuchsbaum.github.io/rMVPA/reference/mock_context.md),
and
[`validate_plugin_model()`](http://bbuchsbaum.github.io/rMVPA/reference/validate_plugin_model.md)
to verify your contract before expensive runs.

``` r
ds <- gen_sample_dataset(c(4, 4, 4), nobs = 20, nlevels = 2, blocks = 2)
mspec <- toy_plugin_model(ds$dataset, ds$design)

plugin_check <- validate_plugin_model(
  mspec,
  roi_data = mock_roi_data(n_train = 20, n_features = 8, n_test = 0, seed = 1),
  context = mock_context(design = ds$design, id = 101L)
)

plugin_check
#> <plugin_validation_result> class=toy_plugin_model, metrics=2
#>   metric_names: score, spread
```

``` r
validate_model_spec(mspec, require_schema = TRUE, dry_run = TRUE)
#> <model_spec_validation_result> class=toy_plugin_model valid=TRUE (pass=5 warn=0 fail=0)
```

## Run on searchlight and regional harnesses

``` r
sl_res <- run_searchlight(mspec, radius = 2, method = "standard")
names(sl_res$results)
#> [1] "score"  "spread"

region_mask <- NeuroVol(
  sample(1:3, size = length(ds$dataset$mask), replace = TRUE),
  space(ds$dataset$mask)
)
reg_res <- run_regional(mspec, region_mask)
names(reg_res$vol_results)
#> [1] "score"  "spread"
```

## Optional: worker-safe object trimming

If your model spec stores heavy transient objects, define a class method
that strips them before distributed execution.

``` r
strip_dataset.toy_plugin_model <- function(obj, ...) {
  obj$dataset <- NULL
  obj
}
```

## Plugin package checklist

1.  Export your constructor.
2.  Provide `fit_roi.<class>` and `output_schema.<class>`.
3.  Keep metric names stable; schema and runtime metrics must match
    exactly.
4.  Add unit tests around
    [`validate_plugin_model()`](http://bbuchsbaum.github.io/rMVPA/reference/validate_plugin_model.md).
5.  Add one lightweight integration test for
    [`run_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
    and/or
    [`run_regional()`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md).
