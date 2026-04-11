# Explain Searchlight Engine Selection

Reports which searchlight engine would be selected for a given model and
method, along with the eligibility status of each registered engine.

## Usage

``` r
explain_searchlight_engine(
  model_spec,
  method = c("standard", "randomized", "resampled"),
  engine = c("auto", "legacy", "swift", "dual_lda_fast")
)
```

## Arguments

- model_spec:

  A model specification.

- method:

  Searchlight method to audit.

- engine:

  Requested engine policy: `"auto"`, `"legacy"`, `"swift"`, or
  `"dual_lda_fast"`.

## Value

A data frame with selection metadata.
