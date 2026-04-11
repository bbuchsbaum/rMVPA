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
if (FALSE) { # \dontrun{
# Create a MANOVA model using gen_sample_dataset
dset <- gen_sample_dataset(D = c(5, 5, 5), nobs = 50, nlevels = 3, blocks = 3)

# Create dissimilarity matrices for MANOVA design
formula <- y ~ x1 + x2
data_list <- list(
  y = matrix(rnorm(9), nrow = 3),
  x1 = matrix(rnorm(9), nrow = 3),
  x2 = matrix(rnorm(9), nrow = 3)
)
design <- manova_design(formula, data_list)
manova_model_obj <- manova_model(dset$dataset, design)
} # }
```
