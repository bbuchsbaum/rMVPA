# Install rMVPA Command-Line Wrappers

Copies the packaged command-line wrappers into a directory on your
`PATH`. The installed commands are `rmvpa-searchlight` and
`rmvpa-regional`.

## Usage

``` r
install_cli(
  dest_dir,
  overwrite = FALSE,
  commands = c("searchlight", "regional")
)
```

## Arguments

- dest_dir:

  Destination directory for the wrappers, for example `"~/.local/bin"`.
  Required: nothing is written to your home directory unless you choose
  it.

- overwrite:

  Logical; overwrite existing wrapper files if `TRUE`.

- commands:

  Which wrappers to install. Any subset of
  `c("searchlight", "regional")`.

## Value

Invisibly, a named character vector of installed wrapper paths.
