# Cross-Validation Availability

Determine whether cross-validation is specified for the object.

## Usage

``` r
has_crossval(obj)
```

## Arguments

- obj:

  Model specification object.

## Value

Logical indicating if cross-validation will be performed.

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10)
cval <- blocked_cross_validation(ds$design$block_var)
mdl <- load_model("sda_notune")
mspec <- mvpa_model(mdl, ds$dataset, ds$design,
                    "classification", crossval = cval)
has_crossval(mspec)
#> [1] TRUE
```
