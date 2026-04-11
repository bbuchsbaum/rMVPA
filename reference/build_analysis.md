# Build an Analysis Object

Converts an `rmvpa_config` into a concrete analysis object with
initialized design, datasets, and model specification(s).

## Usage

``` r
build_analysis(config)
```

## Arguments

- config:

  An `rmvpa_config` returned by
  [`mvpa_config`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_config.md),
  or a named list containing at least `mode`.

## Value

An object of class `c("rmvpa_analysis", "list")`.
