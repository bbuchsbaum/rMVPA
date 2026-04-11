# Training Labels/Response Extraction

Extract the training labels or response variable from an object.

## Usage

``` r
y_train(obj)

# S3 method for class 'mvpa_design'
y_train(obj)

# S3 method for class 'item_design'
y_train(obj)

# S3 method for class 'mvpa_model'
y_train(obj)

# S3 method for class 'model_spec'
y_train(obj)

# S3 method for class 'feature_rsa_model'
y_train(obj)

# S3 method for class 'feature_rsa_design'
y_train(obj)

# S3 method for class 'hrfdecoder_model'
y_train(obj)
```

## Arguments

- obj:

  The object from which to extract the training response variable.

## Value

The training response variable (factor for classification, numeric for
regression).

## Examples

``` r
ds <- gen_sample_dataset(D = c(4, 4, 4), nobs = 10)
y_train(ds$design)
#>  [1] c d a b b d e e c a
#> Levels: a b c d e
```
