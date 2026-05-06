# Extending rMVPA: Plugin Analysis Models

rMVPA can be extended from external packages by defining S3 model
classes that inherit from `model_spec`. The preferred extension path is
the modern `fit_roi`/`output_schema` contract.

## Recommended Contract

1\. \*\*Constructor\*\* Build a model object with
[`create_model_spec`](http://bbuchsbaum.github.io/rMVPA/reference/create_model_spec.md),
assigning your model class name (for example `"my_model"`).

2\. \*\*`fit_roi.my_model` (required)\*\* Implement per-ROI work and
return
[`roi_result`](http://bbuchsbaum.github.io/rMVPA/reference/roi_result.md).

3\. \*\*`output_schema.my_model` (recommended)\*\* Return a named list
with values `"scalar"` or `"vector[N]"`. This enables schema-driven map
construction in searchlight analyses.

4\. \*\*Optional overrides\*\* - `y_train.my_model` / `y_test.my_model`
if labels are not available from the attached design object. -
`strip_dataset.my_model` for advanced parallel memory control.

## Plugin Test Helpers

\-
[`mock_roi_data`](http://bbuchsbaum.github.io/rMVPA/reference/mock_roi_data.md)
builds lightweight ROI payloads for unit tests. -
[`mock_context`](http://bbuchsbaum.github.io/rMVPA/reference/mock_context.md)
builds the standard `fit_roi` context list. -
[`validate_plugin_model`](http://bbuchsbaum.github.io/rMVPA/reference/validate_plugin_model.md)
runs a contract check on one ROI and validates schema/metric
agreement. -
[`validate_model_spec`](http://bbuchsbaum.github.io/rMVPA/reference/validate_model_spec.md)
performs a plugin-readiness lint with optional dry-run execution.

## Two Registry Types

rMVPA exposes two different extension registries:

- Classifier registry:

  Use
  [`register_mvpa_model`](http://bbuchsbaum.github.io/rMVPA/reference/register_mvpa_model.md)
  to add a new low-level estimator (fit/predict/prob) that plugs into
  the built-in `mvpa_model` analysis type.

- Analysis-type registry:

  Define a new S3 model class (via
  [`create_model_spec`](http://bbuchsbaum.github.io/rMVPA/reference/create_model_spec.md) +
  `fit_roi.<class>`) when you need a new per-ROI analysis workflow, not
  just a new estimator backend.

## Notes

\-
[`process_roi.default`](http://bbuchsbaum.github.io/rMVPA/reference/process_roi-methods.md)
automatically prefers `fit_roi` when a method exists for your class. -
Legacy `train_model`/`process_roi` extension paths are still available
for backward compatibility, but new plugins should use `fit_roi`.

## Minimal Example


    my_model <- function(dataset, design, alpha = 0.1, ...) {
      create_model_spec(
        "my_model",
        dataset = dataset,
        design = design,
        alpha = alpha,
        ...
      )
    }

    fit_roi.my_model <- function(model, roi_data, context, ...) {
      score <- mean(roi_data$train_data)
      roi_result(
        metrics = c(score = score),
        indices = roi_data$indices,
        id = context$id
      )
    }

    output_schema.my_model <- function(model) {
      list(score = "scalar")
    }

    # optional only when you do not rely on design defaults:
    # y_train.my_model <- function(obj) y_train(obj$design)
