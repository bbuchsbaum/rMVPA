# Build a Public Analysis Configuration

Creates a lightweight configuration object for the supported high-level
workflow API. Config values can be supplied directly as named arguments
or loaded from a YAML or R config file via `config`.

## Usage

``` r
mvpa_config(mode = c("searchlight", "regional", "global"), config = NULL, ...)
```

## Arguments

- mode:

  Analysis mode: `"searchlight"`, `"regional"`, or `"global"`.

- config:

  Optional path to a YAML or R config file.

- ...:

  Additional named configuration values.

## Value

An object of class `c("rmvpa_config", "list")`.
