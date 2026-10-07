# Summary Method for Feature RSA Model

Summary Method for Feature RSA Model

## Usage

``` r
# S3 method for class 'feature_rsa_model'
summary(object, ...)
```

## Arguments

- object:

  The feature RSA model

- ...:

  Additional args

## Value

A list of summary statistics for the feature RSA model (printed as side
effect).

## Examples

``` r
sample_ds <- gen_sample_dataset(c(4, 4, 4), nobs = 12, blocks = 2)
des <- feature_rsa_design(F = matrix(rnorm(12 * 4), 12, 4),
                          labels = paste0("t", 1:12),
                          block_var = sample_ds$design$block_var)
mdl <- feature_rsa_model(sample_ds$dataset, des, method = "pca")
summary(mdl)
#> ================================================== 
#>           Feature RSA Model           
#> ================================================== 
#> 
#> Method:  pca 
#> Feature standardization:  scale 
#> Return predictions:      No 
#> Number of Observations:  12 
#> Feature Dimensions:      4 
#> Max components limit:    4 
#> Component selection:     loo 
#> Component objective:     mse 
#> One-SE rule:             Yes 
#> Status:  Model not yet trained
#> Cross-Validation:  Configured
#> 
#>  ================================================== 
```
