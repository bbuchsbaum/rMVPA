# Print mvpa_sysinfo Object

Formats and prints the system information gathered by `mvpa_sysinfo`.
This method provides a user-friendly display of the system
configuration.

## Usage

``` r
# S3 method for class 'mvpa_sysinfo'
print(x, ...)
```

## Arguments

- x:

  An object of class \`mvpa_sysinfo\`.

- ...:

  Ignored.

## Value

Invisibly returns the input object `x` (called for side effects).

## Examples

``` r
if (FALSE) { # \dontrun{
  info <- mvpa_sysinfo()
  print(info)
} # }
```
