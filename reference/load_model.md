# Load a Pre-defined MVPA Model

Retrieves a model specification from the pre-defined set of MVPA models.

## Usage

``` r
load_model(name)
```

## Arguments

- name:

  Character string specifying the model to load. Must be a pre-defined
  MVPA model name:

  corclass

  :   Correlation-based classifier with template matching

  sda_notune

  :   Simple Shrinkage Discriminant Analysis without tuning

  sda_boot

  :   SDA with bootstrap resampling

  glmnet_opt

  :   Elastic net with EPSGO parameter optimization

  sparse_sda

  :   SDA with sparsity constraints

  sda_ranking

  :   SDA with automatic feature ranking

  mgsda

  :   Multi-Group Sparse Discriminant Analysis

  lda_thomaz

  :   Modified LDA for high-dimensional data

  hdrda

  :   High-Dimensional Regularized Discriminant Analysis

  spacenet_tvl1

  :   Spatially-regularized sparse linear model with TV-L1 prior
      (SpaceNet-style)

## Value

A list containing the model specification with the following components:

- type:

  Model type: "Classification" or "Regression"

- library:

  Required R package(s) for the model

- label:

  Human-readable model name

- parameters:

  Data frame describing tunable parameters

- grid:

  Function to generate parameter tuning grid

- fit:

  Function to fit the model

- predict:

  Function to generate predictions

- prob:

  Function to generate class probabilities (classification only)

## References

Gramfort, A., Thirion, B., & Varoquaux, G. (2013). *Identifying
predictive regions from fMRI with TV-L1 prior*. Pattern Recognition in
Neuroimaging (PRNI), IEEE. https://inria.hal.science/hal-00839984

## See also

[`MVPAModels`](http://bbuchsbaum.github.io/rMVPA/reference/MVPAModels.md)
for the complete list of available custom MVPA models

[`mvpa_model`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_model.md)
for using these models in MVPA analyses

## Examples

``` r
# Load custom MVPA model
model <- load_model("sda_notune")

# Load correlation classifier with parameter tuning options
corr_model <- load_model("corclass")
print(corr_model$parameters)  # View tunable parameters
#>   parameters     class                                           label
#> 1     method character correlation type: pearson, spearman, or kendall
#> 2     robust   logical                                   mean or huber
```
