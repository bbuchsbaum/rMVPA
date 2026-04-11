# Default method for run_searchlight

By default, if an object's class does not implement a specific
`run_searchlight.<class>` method, this fallback will call
`run_searchlight_base`.

## Usage

``` r
# Default S3 method
run_searchlight(
  model_spec,
  radius = 8,
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

  Additional arguments passed through:

  engine

  :   Searchlight engine: `"auto"` (default), `"legacy"`, `"swift"`, or
      `"dual_lda_fast"`.

  combiner

  :   Combiner function or name for randomized/resampled methods
      (default `"average"`).

  drop_probs

  :   If `TRUE`, drop probability predictions (default `FALSE`).

  fail_fast

  :   If `TRUE`, stop on first ROI error (default `FALSE`).

  k

  :   Number of neighbors for searchlight construction (default `NULL`).

  verbose

  :   If `TRUE`, print progress messages (default `FALSE`).

  incremental, gamma

  :   Engine-specific arguments for dual LDA engine.

## Value

A \`searchlight_result\` object containing spatial maps for each metric.

## Examples

``` r
if (FALSE) { # \dontrun{
  # See run_searchlight generic for examples
} # }
```
