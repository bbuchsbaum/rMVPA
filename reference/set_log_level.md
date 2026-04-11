# Set rMVPA Logging Level

Convenience helper to control the verbosity of rMVPA's logging output.
Internally this wraps
[`futile.logger::flog.threshold()`](https://rdrr.io/pkg/futile.logger/man/flog.threshold.html)
for the package's logger.

## Usage

``` r
set_log_level(level = "INFO")
```

## Arguments

- level:

  Logging level to use. Can be a character string (`"TRACE"`, `"DEBUG"`,
  `"INFO"`, `"WARN"`, `"ERROR"`, `"FATAL"`, `"OFF"`) or a numeric level
  constant from the futile.logger package.

## Value

Invisibly returns the numeric log level.

## Details

Typical usage:

- `set_log_level("INFO")` - default, hides debug messages.

- `set_log_level("DEBUG")` - show detailed per-batch timing and ROI
  diagnostics.

- `set_log_level("WARN")` - only warnings and errors.

This affects all rMVPA logging performed via futile.logger.

## Examples

``` r
if (FALSE) { # \dontrun{
  rMVPA::set_log_level("DEBUG")
  rMVPA::set_log_level("WARN")
} # }
```
