# A "base" function for searchlight analysis

This function implements the generic logic for running a searchlight:

1.  Checks `radius` and `method`.

2.  For "standard" searchlight, calls `do_standard(...)`.

3.  For "randomized", calls `do_randomized(...)` with `niter` times.

4.  Handles the `combiner` function or string ("pool", "average").

## Usage

``` r
run_searchlight_base(
  model_spec,
  radius = 8,
  method = c("standard", "randomized", "resampled"),
  niter = 4,
  combiner = "average",
  drop_probs = FALSE,
  fail_fast = FALSE,
  preflight = c("warn", "error", "off"),
  engine = NULL,
  k = NULL,
  backend = c("default", "shard", "auto"),
  verbose = FALSE,
  ...
)
```

## Arguments

- model_spec:

  A model specification object (e.g., `mvpa_model`, `vector_rsa_model`,
  etc.).

- radius:

  Numeric searchlight radius (1 to 100).

- method:

  Character: "standard" or "randomized".

- niter:

  Number of iterations if `method="randomized"`.

- combiner:

  Either a function that combines partial results or a string ("pool",
  "average") that selects a built-in combiner.

- drop_probs:

  Logical; if TRUE, drop per-ROI probability matrices after computing
  metrics to save memory. Default FALSE.

- fail_fast:

  Logical; if TRUE, stop immediately on first ROI error instead of
  continuing. Default FALSE.

- k:

  Optional integer; number of cross-validation folds. Default NULL (use
  model default).

- ...:

  Additional arguments passed on to `do_standard` or `do_randomized`.

## Value

The result object from `do_standard` or `do_randomized` (often a
`searchlight_result` or similar).

## Details

It does not assume any specific model type, but expects that
`model_spec` is compatible with `do_standard(...)` or
`do_randomized(...)` in your code.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Internal base function - users should call run_searchlight instead
  # result <- run_searchlight_base(model_spec, radius=8, method="standard")
} # }
```
