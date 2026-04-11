# Cross-validation samples

Apply a cross-validation scheme to split the data into training and
testing sets.

## Usage

``` r
crossval_samples(obj, data, y, ...)

# S3 method for class 'sequential_blocked_cross_validation'
crossval_samples(obj, data, y, ...)

# S3 method for class 'kfold_cross_validation'
crossval_samples(obj, data, y, ...)

# S3 method for class 'blocked_cross_validation'
crossval_samples(obj, data, y, ...)

# S3 method for class 'bootstrap_blocked_cross_validation'
crossval_samples(obj, data, y, ...)

# S3 method for class 'custom_cross_validation'
crossval_samples(obj, data, y, id = ".id", ...)

# S3 method for class 'twofold_blocked_cross_validation'
crossval_samples(obj, data, y, ...)

# S3 method for class 'mvpa_model'
crossval_samples(obj, data, y, ...)
```

## Arguments

- obj:

  A cross-validation control object.

- data:

  A data frame containing the predictors.

- y:

  A vector containing the response variable.

- ...:

  Extra arguments passed to the specific cross-validation methods (e.g.,
  \`id\` for custom cross-validation).

- id:

  Column name used for the fold identifier column in the returned
  tibble.

## Value

A tibble containing training and testing sets for each fold.

## Examples

``` r
cval <- kfold_cross_validation(len = 20, nfolds = 4)
dat  <- as.data.frame(matrix(rnorm(20 * 2), 20, 2))
y    <- factor(rep(letters[1:4], 5))
crossval_samples(cval, dat, y)
#> # A tibble: 4 × 5
#>   ytrain       ytest        train               test               .id  
#>   <named list> <named list> <named list>        <named list>       <chr>
#> 1 <fct [15]>   <fct [5]>    <resample [15 x 2]> <resample [5 x 2]> 01   
#> 2 <fct [15]>   <fct [5]>    <resample [15 x 2]> <resample [5 x 2]> 02   
#> 3 <fct [15]>   <fct [5]>    <resample [15 x 2]> <resample [5 x 2]> 03   
#> 4 <fct [15]>   <fct [5]>    <resample [15 x 2]> <resample [5 x 2]> 04   
```
