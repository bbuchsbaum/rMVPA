# Run Searchlight Analysis

Execute a searchlight analysis using multivariate pattern analysis.

## Usage

``` r
run_searchlight(
  model_spec,
  radius,
  method = c("standard", "randomized", "resampled"),
  niter = 4,
  backend = c("default", "shard", "auto"),
  ...
)

# S3 method for class 'feature_rsa_model'
run_searchlight(
  model_spec,
  radius,
  method = c("standard", "randomized", "resampled"),
  niter = 4,
  backend = c("default", "shard", "auto"),
  ...
)
```

## Arguments

- model_spec:

  A `mvpa_model` instance containing the model specifications

- radius:

  The searchlight radius in millimeters

- method:

  The type of searchlight, either 'randomized' or 'standard'

- niter:

  The number of searchlight iterations (used only for 'randomized'
  method)

- backend:

  Execution backend: `"default"` (standard pipeline), `"shard"`
  (shared-memory backend), or `"auto"` (try shard and fall back to
  default).

- ...:

  Extra arguments passed to specific searchlight methods. Currently
  supported:

  - `batch_size`: Integer specifying how many searchlights (ROIs) are
    grouped into a single batch for processing by
    [`mvpa_iterate`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md).
    Each batch is launched sequentially, while ROIs within a batch are
    processed in parallel (using the active future/furrr plan). The
    default is 10% of the total number of searchlights. Smaller values
    lower peak memory usage and can improve load balancing, at the cost
    of more scheduling overhead; larger values reduce overhead but
    require more memory per worker. As a rule of thumb, start with the
    default, decrease `batch_size` if you hit memory limits, and
    increase it if you have many ROIs, ample RAM, and see low CPU
    utilization.

## Value

A named list of `NeuroVol` objects containing performance metrics (e.g.,
AUC) at each voxel location

## Progress reporting

When `verbose = TRUE` is passed via `...` and the progressr package is
installed, real-time per-ROI progress updates are emitted from parallel
workers. To enable a progress bar:


      library(progressr)
      handlers(global = TRUE)               # once per session
      result <- run_searchlight(mspec, radius = 8, verbose = TRUE)

Without progressr, only coarse batch-level log messages are shown.

## Examples

``` r
# \donttest{
  # Generate sample dataset with categorical response
  dataset <- gen_sample_dataset(
    D = c(5,5,5),           # 5x5x5 volume
    nobs = 60,              # 60 observations
    response_type = "categorical",
    data_mode = "image",
    blocks = 3,             # 3 blocks for cross-validation
    nlevels = 2             # binary classification
  )
  
  # Create cross-validation specification using blocks
  cval <- blocked_cross_validation(dataset$design$block_var)
  
  # Load the SDA classifier (Shrinkage Discriminant Analysis)
  model <- load_model("sda_notune")
  
  # Create MVPA model
  mspec <- mvpa_model(
    model = model,
    dataset = dataset$dataset,
    design = dataset$design,
    model_type = "classification",
    crossval = cval
  )
  
  # Run searchlight analysis
  results <- run_searchlight(
    mspec,
    radius = 3,            # 3mm radius
    method = "standard"    # Use standard searchlight
  )
#> INFO [2026-10-08 02:06:39] searchlight engine: sda_fast
  
  # Run with custom batch size for memory management
  # results <- run_searchlight(
  #   mspec,
  #   radius = 3,
  #   method = "standard",
  #   batch_size = 500      # Process 500 searchlights per batch
  # )
  
  # Results contain performance metrics
  # Access them with results$performance
# }
```
