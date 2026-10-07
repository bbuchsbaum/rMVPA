# Create a MANOVA Model

This function creates a MANOVA model object containing an
\`mvpa_dataset\` instance and a \`manova_design\` instance.

## Usage

``` r
manova_model(dataset, design)
```

## Arguments

- dataset:

  An `mvpa_dataset` instance.

- design:

  A `manova_design` instance.

## Value

A MANOVA model object with class attributes "manova_model" and "list".

## Details

The function takes an \`mvpa_dataset\` instance and a \`manova_design\`
instance as input, and returns a MANOVA model object. The object is a
list that contains the dataset and the design with class attributes
"manova_model" and "list". This object can be used for further
multivariate statistical analysis using the MANOVA method.

## Examples

``` r
# Create a MANOVA model using gen_sample_dataset
dset <- gen_sample_dataset(D = c(4, 4, 4), nobs = 24, nlevels = 3, blocks = 3)

# MANOVA design: voxel patterns modelled by condition and block
design <- manova_design(~ Y + block_var, dset$design$train_design)
manova_model_obj <- manova_model(dset$dataset, design)
```
