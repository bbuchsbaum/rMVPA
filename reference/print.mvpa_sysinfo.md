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
# mvpa_sysinfo() prints on creation; capture that output, then print explicitly
invisible(utils::capture.output(info <- mvpa_sysinfo()))
print(info)
#> ------------------------- rMVPA System Information -------------------------
#> R Version                : N/A
#> Platform                 : x86_64-pc-linux-gnu
#> Operating System         : Linux 6.17.0-1022-azure
#> Node Name                : runnervmmprz5
#> User                     : runner
#> Locale                   : LC_CTYPE=C.UTF-8;LC_NUMERIC=C;LC_TIME=C.UTF-8;LC_COLLATE=C;LC_MONETARY=C.UTF-8;LC_MESSAGES=C.UTF-8;LC_PAPER=C.UTF-8;LC_NAME=C;LC_ADDRESS=C;LC_TELEPHONE=C;LC_MEASUREMENT=C.UTF-8;LC_IDENTIFICATION=C
#> rMVPA Version            : 0.1.3
#> Parallel Backend (Future): sequential
#> 
#> Key Dependencies:
#>   - neuroim2 : 0.19.1
#>   - neurosurf: 0.1.0
#>   - rsample  : 1.3.2
#>   - yardstick: 1.4.0
#>   - future   : 1.76.0
#>   - furrr    : 0.4.0
#>   - dplyr    : 1.2.1
#>   - tibble   : 3.3.1
#>   - purrr    : 1.2.2
#>   - stats    : 4.6.1
#>   - MASS     : 7.3.65
#> --------------------------------------------------------------------------
#> 
```
