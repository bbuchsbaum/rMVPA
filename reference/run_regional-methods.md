# Region of Interest Based MVPA Analysis

Run a separate MVPA analysis for multiple disjoint regions of interest.

## Usage

``` r
run_regional(
  model_spec,
  region_mask,
  backend = c("default", "shard", "auto"),
  ...
)

run_regional_base(
  model_spec,
  region_mask,
  coalesce_design_vars = FALSE,
  processor = NULL,
  verbose = FALSE,
  compute_performance = model_spec$compute_performance,
  return_predictions = model_spec$return_predictions,
  return_fits = model_spec$return_fits,
  pool_predictions = c("none", "mean", "stack"),
  pooled_weights = NULL,
  stack_folds = NULL,
  stack_seed = NULL,
  stack_lambda = 0.001,
  save_rdm_vectors_dir = NULL,
  preflight = c("warn", "error", "off"),
  backend = c("default", "shard", "auto"),
  ...
)

# Default S3 method
run_regional(
  model_spec,
  region_mask,
  backend = c("default", "shard", "auto"),
  ...
)

# S3 method for class 'mvpa_model'
run_regional(
  model_spec,
  region_mask,
  backend = c("default", "shard", "auto"),
  ...
)

# S3 method for class 'rsa_model'
run_regional(
  model_spec,
  region_mask,
  backend = c("default", "shard", "auto"),
  ...
)

# S3 method for class 'vector_rsa_model'
run_regional(
  model_spec,
  region_mask,
  backend = c("default", "shard", "auto"),
  ...
)
```

## Arguments

- model_spec:

  A `mvpa_model` instance containing the model specifications

- region_mask:

  A `NeuroVol` or `NeuroSurface` object where each region is identified
  by a unique integer

- backend:

  Execution backend: `"default"` (standard pipeline), `"shard"`
  (shared-memory backend), or `"auto"` (try shard and fall back to
  default).

- ...:

  Extra arguments passed to specific regional analysis methods (e.g.,
  \`return_fits\`, \`compute_performance\`).

- coalesce_design_vars:

  If `TRUE`, merges design variables into the prediction table (if
  present and generated). Default is `FALSE`.

- processor:

  An optional custom processor function for each region (ROI). If NULL
  (default), behavior depends on the `model_spec` class.

- verbose:

  If `TRUE`, print progress messages during iteration (default is
  `FALSE`).

- compute_performance:

  Logical indicating whether to compute performance metrics (default
  `TRUE`).

- return_predictions:

  Logical indicating whether to combine a full prediction table
  (defaults to `model_spec$return_predictions`).

- return_fits:

  Logical indicating whether to return the fitted models (default
  `FALSE`).

- pool_predictions:

  Character scalar controlling pooled outputs: `"none"` (default),
  `"mean"` (weighted mean pooling across ROIs), or `"stack"` (cheap
  cross-fitted stacking on OOF ROI predictions).

- pooled_weights:

  Optional numeric vector of ROI weights used when
  `pool_predictions = "mean"`. Must have one weight per ROI.

- stack_folds:

  Optional fold specification for stacking. Can be: an integer number of
  folds, a fold-id vector, or a list of fold indices.

- stack_seed:

  Optional seed used when auto-generating stacking folds.

- stack_lambda:

  Ridge penalty used by glmnet in stacking.

- save_rdm_vectors_dir:

  Optional directory where feature-RSA predicted and observed RDM
  vectors should be written batch-by-batch. When supplied for
  \`feature_rsa_model(..., return_rdm_vectors = TRUE)\`,
  \`run_regional()\` writes compact \`rdm_batches/batch\_\*.rds\` files
  in that directory and returns \`rdm_batch_dir\` in the result instead
  of retaining all ROI RDM vectors in memory.

- preflight:

  One of `"warn"` (default), `"error"`, or `"off"` controlling whether
  analysis preflight findings emit warnings, stop the run, or are
  skipped.

## Value

A `regional_mvpa_result` object (list) containing:

- performance_table:

  A tibble of performance metrics for each region (if computed).

- prediction_table:

  A tibble with detailed predictions for each observation/region (if
  generated).

- pooled_prediction_table:

  Optional pooled trial-level prediction table when pooling is
  requested.

- pooled_performance:

  Optional pooled performance metrics when pooling is requested.

- vol_results:

  A list of volumetric maps representing performance metrics across
  space (if computed).

- fits:

  A list of fitted model objects for each region (if requested via
  \`return_fits=TRUE\`).

- rdm_batch_dir:

  Optional directory of file-backed feature-RSA RDM-vector batches when
  \`save_rdm_vectors_dir\` is used.

- model_spec:

  The original model specification object provided.

\# Note: Original documentation said 'performance', clarified here.

## Details

This function serves as the base implementation for regional analyses,
orchestrating data preparation, iteration over regions, performance
computation, and result aggregation. Specific \`run_regional\` methods
for different model classes may call this function or provide
specialized behavior.

This is the fallback method called when no specialized \`run_regional\`
method is found for the class of \`model_spec\`. It typically calls
\`run_regional_base\`.

This method provides the standard regional analysis pipeline for objects
of class \`mvpa_model\` by calling \`run_regional_base\`.

For \`rsa_model\` objects, \`return_predictions\` defaults to \`FALSE\`
as standard RSA typically doesn't produce a prediction table in the same
way as classification/regression models.

For \`vector_rsa_model\` objects, \`return_predictions\` defaults to
\`FALSE\` in \`run_regional_base\`. If
\`model_spec\$return_predictions\` is TRUE, this method will assemble an
\`observation_scores_table\`.

## Progress reporting

When `verbose = TRUE` is passed via `...` and the progressr package is
installed, real-time per-ROI progress updates are emitted from parallel
workers. To enable a progress bar:


      library(progressr)
      handlers(global = TRUE)               # once per session
      result <- run_regional(mspec, region_mask, verbose = TRUE)

Without progressr, only coarse batch-level log messages are shown.

## Examples

``` r
# \donttest{
  # Generate sample dataset (3D volume with categorical response)
  dataset <- gen_sample_dataset(
    D = c(10,10,10),       # Small 10x10x10 volume
    nobs = 100,            # 100 observations
    nlevels = 3,           # 3 classes
    response_type = "categorical",
    data_mode = "image",
    blocks = 3             # 3 blocks for cross-validation
  )
  
  # Create region mask with 5 ROIs
  region_mask <- neuroim2::NeuroVol(
    sample(1:5, size=length(dataset$dataset$mask), replace=TRUE),
    neuroim2::space(dataset$dataset$mask)
  )
  
  # Create cross-validation specification
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  # Load SDA classifier (Shrinkage Discriminant Analysis)
  model <- load_model("sda_notune")
  
  # Create MVPA model
  mspec <- mvpa_model(
    model = model,
    dataset = dataset$dataset,
    design = dataset$design,
    model_type = "classification",
    crossval = cval,
    return_fits = TRUE    # Return fitted models
  )
  
  # Run regional analysis
  results <- run_regional(mspec, region_mask)
#> INFO [2026-05-08 14:27:43] 
#> MVPA Iteration Complete
#> - Total ROIs: 5
#> - Processed: 5
#> - Skipped: 0
#> INFO [2026-05-08 14:27:43] run_regional: 5 ROIs processed (success=5, errors=0)
  
  # Access results
  head(results$performance_table)     # Performance metrics
#> # A tibble: 5 × 3
#>   roinum Accuracy     AUC
#>    <int>    <dbl>   <dbl>
#> 1      1     0.32 -0.0887
#> 2      2     0.39  0.0772
#> 3      3     0.36  0.0554
#> 4      4     0.31 -0.0872
#> 5      5     0.31 -0.105 
  head(results$prediction_table)      # Predictions
#> # A tibble: 6 × 9
#> # Rowwise: 
#>   .rownum roinum observed pobserved predicted correct   prob_a prob_b  prob_c
#>     <int>  <int> <fct>        <dbl> <chr>     <lgl>      <dbl>  <dbl>   <dbl>
#> 1       1      1 c         0.169    b         FALSE   0.0458    0.785 0.169  
#> 2       2      1 c         0.0669   a         FALSE   0.482     0.451 0.0669 
#> 3       3      1 a         0.000284 b         FALSE   0.000284  0.995 0.00455
#> 4       4      1 b         0.168    a         FALSE   0.790     0.168 0.0422 
#> 5       5      1 c         0.0639   a         FALSE   0.507     0.429 0.0639 
#> 6       6      1 b         0.442    b         TRUE    0.245     0.442 0.313  
  first_roi_fit <- results$fits[[1]]  # First ROI's fitted model
# }
```
