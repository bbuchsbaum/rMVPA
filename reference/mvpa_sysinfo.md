# Report System and Package Information for rMVPA

Gathers and displays information about the R session, operating system,
rMVPA version, and key dependencies. This information is helpful for
debugging, reporting issues, and ensuring reproducibility.

## Usage

``` r
mvpa_sysinfo()
```

## Value

Invisibly returns a list containing the gathered system and package
information. It is primarily called for its side effect: printing the
formatted information to the console.

## Examples

``` r
# Display system information and capture it in a variable
sys_info <- mvpa_sysinfo()
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
sys_info$platform
#> [1] "x86_64-pc-linux-gnu"
sys_info$dependencies$rsample
#> [1] "1.3.2"
```
