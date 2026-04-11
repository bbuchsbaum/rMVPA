# Pre-defined MVPA Classification Models

An environment containing custom classification models for MVPA
analysis.

## Usage

``` r
MVPAModels
```

## Format

An environment with the following models:

- corclass:

  Correlation-based classifier using template matching with options
  (pearson, spearman, kendall)

- corsim:

  Alias for corclass

- sda_notune:

  Shrinkage Discriminant Analysis (SDA) without parameter tuning

- sda_boot:

  SDA with bootstrap resampling and feature selection

- glmnet_opt:

  Elastic net classifier (glmnet) with optimized alpha/lambda via EPSGO

- sparse_sda:

  SDA with sparsity constraints and feature selection

- sda_ranking:

  SDA with feature ranking and selection via higher criticism

- mgsda:

  Multi-Group Sparse Discriminant Analysis

- lda_thomaz:

  Modified LDA classifier for high-dimensional data

- hdrda:

  High-Dimensional Regularized Discriminant Analysis

- spacenet_tvl1:

  Spatially-regularized sparse linear model with TV-L1 penalty for
  global whole-brain analysis

## Value

An environment containing registered MVPA model specifications.

## Details

Models are accessed via `load_model(name)`. Each model specification
includes `fit`, `predict`, and `prob` methods.

The `spacenet_tvl1` model follows the SpaceNet formulation used in
Nilearn.

## References

Gramfort, A., Thirion, B., & Varoquaux, G. (2013). *Identifying
predictive regions from fMRI with TV-L1 prior*. Pattern Recognition in
Neuroimaging (PRNI), IEEE. https://inria.hal.science/hal-00839984

## Examples

``` r
# Load simple SDA classifier
model <- load_model("sda_notune")

# Load correlation classifier
model <- load_model("corclass")
```
