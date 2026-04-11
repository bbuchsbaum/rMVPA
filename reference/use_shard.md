# Enable Shared-Memory Backend for an MVPA Model Spec

Prepares the dataset for shared-memory access and tags the model
specification so that
[`mvpa_iterate`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md)
uses the shard backend instead of the default furrr pipeline.

## Usage

``` r
use_shard(mod_spec)
```

## Arguments

- mod_spec:

  A model specification created by
  [`mvpa_model`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md)
  (or any constructor that produces an object inheriting from
  `"model_spec"`).

## Value

A copy of `mod_spec` with class `"shard_model_spec"` prepended and a
`$shard_data` list attached.

## Details

This is an **experimental** feature. It requires the shard package and a
platform that supports POSIX shared memory (`shm_open`).

## Examples

``` r
if (FALSE) { # \dontrun{
  mspec <- mvpa_model(mdl, dataset, design, "classification", crossval = cval)
  mspec <- use_shard(mspec)
  results <- run_searchlight(mspec, ...)
} # }
```
