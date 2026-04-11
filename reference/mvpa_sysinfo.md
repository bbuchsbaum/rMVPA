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
if (FALSE) { # \dontrun{
# Display system information in the console
mvpa_sysinfo()

# Capture the information in a variable
sys_info <- mvpa_sysinfo()
print(sys_info$r_version)
print(sys_info$dependencies$rsample)
} # }
```
