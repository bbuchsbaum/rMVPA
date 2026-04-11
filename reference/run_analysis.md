# Run a Built Analysis

Runs a high-level workflow object created by
[`build_analysis`](http://bbuchsbaum.github.io/rMVPA/reference/build_analysis.md).
The input can also be an `rmvpa_config`, in which case it is built
first.

## Usage

``` r
run_analysis(x, preflight = c("warn", "error", "off"), ...)
```

## Arguments

- x:

  An `rmvpa_analysis` or `rmvpa_config`.

- preflight:

  One of `"warn"` (default), `"error"`, or `"off"`.

- ...:

  Additional arguments forwarded to the underlying runner.

## Value

An object of class `c("rmvpa_analysis_run", "list")`.
