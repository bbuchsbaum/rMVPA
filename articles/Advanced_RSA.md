# Advanced RSA Methods: Feature-Based and Vector-Based Approaches

## Introduction

Standard Representational Similarity Analysis (RSA) compares neural
response patterns with model‑based similarity structures. In practice,
you may want stronger control over the feature space, or more efficient
handling of distance vectors and block structure. This vignette
introduces two specialized variants in `rMVPA`—Feature‑Based RSA and
Vector‑Based RSA—and shows when they offer advantages over standard RSA.

## Feature-Based RSA

### Overview

Feature‑Based RSA is useful when you have a well‑defined feature space
for your stimuli and want to directly map neural activity into that
space—for example, reconstructing semantic or visual features from brain
responses. Rather than correlating similarity matrices, the method
learns a mapping from neural patterns to feature dimensions and
evaluates how well those features are recovered.

### Key differences from standard RSA

The core distinction is that Feature‑Based RSA predicts feature values
directly, whereas standard RSA compares similarity structures. Outputs
are therefore different: standard RSA produces associations between
RDMs, while Feature‑Based RSA yields predicted feature values per
stimulus. This makes Feature‑Based RSA a natural choice for
reconstruction tasks, while standard RSA remains ideal for testing
theoretical representational models.

### Implementation Example

Let’s walk through a complete example:

``` r

## Generate data and weakly-correlated features
# Synthetic dataset: 6x6x6, 50 observations (4 blocks)
data_out <- rMVPA::gen_sample_dataset(D = c(6,6,6), nobs = 50, blocks = 4, nlevels = 2)
print(data_out)
#> $dataset
#> 
#>  MVPA Dataset 
#> 
#> - Training Data 
#>   - Dimensions:  6 x 6 x 6 x 50 observations 
#>   - Type:  DenseNeuroVec 
#> - Test Data 
#>   -  None 
#> - Mask Information 
#>   - Areas:  TRUE : 216 
#>   - Active voxels/vertices:  216 
#> 
#> 
#> $design
#> 
#>  MVPA Design 
#> 
#> - Training Data 
#>   - Observations:  50 
#>   - Response Type:  Factor
#>   - Levels:  a, b 
#>   - Class Distribution:  a: 25, b: 25 
#> - Test Data 
#>   -  None 
#> - Structure 
#>   - Blocking:  Present
#>   - Number of Blocks:  4 
#>   - Mean Block Size:  12  (SD:  0.58 ) 
#>   - Split Groups:  None

set.seed(123)
n_stimuli  <- 50
n_features <- 5

# Latent factors shared between features and brain data (weak signal)
Z <- matrix(rnorm(n_stimuli * 2), n_stimuli, 2)

# Build a feature matrix with moderate correlation: F = 0.7 * (Z %*% B) + 0.3 * noise
B <- matrix(rnorm(2 * n_features), 2, n_features)
feature_matrix <- 0.7 * (Z %*% B) + 0.3 * matrix(rnorm(n_stimuli * n_features), n_stimuli, n_features)
feature_matrix <- base::scale(feature_matrix)
colnames(feature_matrix) <- paste0("feature_", 1:n_features)

# Create stimulus labels
stim_labels <- paste0("stim_", 1:n_stimuli)

# Create a subtle latent signal in the brain data driven by the same Z
mask_vol <- data_out$dataset$mask
# Extract logical mask values and indices
mask_vals <- as.logical(neuroim2::values(mask_vol))
mask_idx <- which(mask_vals)

# Extract masked voxel-by-time matrix from NeuroVec (V x T)
vol_list <- neuroim2::vols(data_out$dataset$train_data)
datamat  <- do.call(cbind, lapply(vol_list, function(v) as.numeric(v[mask_idx])))

# Two random spatial patterns within the mask
set.seed(124)
p1 <- rnorm(length(mask_idx))
p2 <- rnorm(length(mask_idx))
p1 <- p1 / sd(p1); p2 <- p2 / sd(p2)

# Inject moderate signal into brain data: X <- X + 0.5 * (p1 %*% t(Z1)) + 0.5 * (p2 %*% t(Z2))
datamat <- datamat + p1 %*% t(0.5 * Z[,1]) + p2 %*% t(0.5 * Z[,2])

# Repackage as a NeuroVec using the existing space and logical mask
train_vec <- neuroim2::SparseNeuroVec(datamat, neuroim2::space(data_out$dataset$train_data), mask = mask_vals)
dset <- mvpa_dataset(train_vec, mask = mask_vol)

# Create feature RSA design from the correlated feature matrix.
# Limit max_comps so PLS and PCA use a low‑dimensional subspace,
# making their behavior distinguishable in the example.
feature_design <- feature_rsa_design(
  F = feature_matrix,
  labels = stim_labels,
  max_comps = 2
)

# Create cross-validation structure using the block information
crossval <- blocked_cross_validation(data_out$design$block_var)

# Create feature RSA model
feature_model <- feature_rsa_model(
  dataset = dset,
  design = feature_design,
  method = "pls",  # Partial Least Squares
  crossval = crossval  # Add cross-validation
)

# Create proper region mask from the dataset's mask
mask_vol <- data_out$dataset$mask
nvox <- sum(mask_vol)
region_mask <- neuroim2::NeuroVol(
  sample(1:3, size = nvox, replace = TRUE),  # 3 regions
  space(mask_vol),
  indices = which(mask_vol > 0)
)

# Run regional analysis (suppress log messages for clarity)
results <- suppressWarnings(suppressMessages(run_regional(feature_model, region_mask)))

# Examine results
print(results$performance_table)
#> # A tibble: 3 × 10
#>   roinum pattern_correlation pattern_discrimination pattern_rank_percentile
#>    <int>               <dbl>                  <dbl>                   <dbl>
#> 1      1               0.403                  0.378                   0.798
#> 2      2               0.385                  0.349                   0.786
#> 3      3               0.368                  0.335                   0.766
#> # ℹ 6 more variables: rdm_correlation <dbl>, voxel_correlation <dbl>,
#> #   mse <dbl>, r_squared <dbl>, mean_voxelwise_temporal_cor <dbl>, ncomp <dbl>
```

### Available methods

Feature‑Based RSA supports several estimators. PLS finds components that
maximize covariance between brain data and features, and works well when
features are correlated. PCA reduces the feature space before
regression, providing a compact, interpretable basis. GLMNet (elastic
net) adds regularization that helps with feature selection and
multicollinearity.

``` r

# Compare different methods
methods <- c("pls", "pca", "glmnet")
results_list <- lapply(methods, function(method) {
  model <- feature_rsa_model(
    dataset = dset,
    design = feature_design,
    method = method,
    crossval = crossval  # Add cross-validation
  )
  suppressWarnings(suppressMessages(run_regional(model, region_mask)))
})

# Compare performance
for (i in seq_along(methods)) {
  cat("\nMethod:", methods[i], "\n")
  print(results_list[[i]]$performance_table)
}
#> 
#> Method: pls 
#> # A tibble: 3 × 10
#>   roinum pattern_correlation pattern_discrimination pattern_rank_percentile
#>    <int>               <dbl>                  <dbl>                   <dbl>
#> 1      1               0.403                  0.378                   0.798
#> 2      2               0.385                  0.349                   0.786
#> 3      3               0.368                  0.335                   0.766
#> # ℹ 6 more variables: rdm_correlation <dbl>, voxel_correlation <dbl>,
#> #   mse <dbl>, r_squared <dbl>, mean_voxelwise_temporal_cor <dbl>, ncomp <dbl>
#> 
#> Method: pca 
#> # A tibble: 3 × 10
#>   roinum pattern_correlation pattern_discrimination pattern_rank_percentile
#>    <int>               <dbl>                  <dbl>                   <dbl>
#> 1      1               0.399                  0.374                   0.794
#> 2      2               0.380                  0.344                   0.779
#> 3      3               0.363                  0.330                   0.760
#> # ℹ 6 more variables: rdm_correlation <dbl>, voxel_correlation <dbl>,
#> #   mse <dbl>, r_squared <dbl>, mean_voxelwise_temporal_cor <dbl>, ncomp <dbl>
#> 
#> Method: glmnet 
#> # A tibble: 3 × 10
#>   roinum pattern_correlation pattern_discrimination pattern_rank_percentile
#>    <int>               <dbl>                  <dbl>                   <dbl>
#> 1      1               0.379                  0.357                   0.803
#> 2      2               0.367                  0.339                   0.787
#> 3      3               0.377                  0.351                   0.793
#> # ℹ 6 more variables: rdm_correlation <dbl>, voxel_correlation <dbl>,
#> #   mse <dbl>, r_squared <dbl>, mean_voxelwise_temporal_cor <dbl>, ncomp <dbl>
```

## Vector-Based RSA

### Overview

Vector‑Based RSA starts from pre‑computed distance matrices and operates
directly on vectorized distances. It is particularly helpful when you
need to ignore within‑block comparisons or scale to large datasets—it
avoids materializing full similarity matrices and handles blocks
efficiently.

### Key differences from standard RSA

Compared to standard RSA, the vector approach (i) works on distance
vectors instead of full matrices, (ii) builds block‑wise exclusions into
the workflow, and (iii) reduces memory use by retaining only the
necessary comparisons. These properties make it a good fit for
high‑resolution or long‑duration datasets.

### Implementation Example

``` r

# Create distance matrix for stimuli
stim_distances <- as.matrix(dist(feature_matrix))
rownames(stim_distances) <- stim_labels

# Create block structure (e.g., runs)
blocks <- rep(1:5, length.out = n_stimuli)

# Create vector RSA design
vector_design <- vector_rsa_design(
    D = stim_distances,
    labels = stim_labels,
    block_var = blocks
)

# Create vector RSA model
vector_model <- vector_rsa_model(
  dataset = dset,
  design = vector_design,
  distfun = cordist(),  # Correlation distance
  rsa_simfun = "pearson"
)

# Run analysis (suppress log messages for clarity)
results_vector <- suppressWarnings(suppressMessages(run_regional(vector_model, region_mask)))

# Examine results
print(results_vector$performance_table)
#> # A tibble: 3 × 2
#>   roinum rsa_score
#>    <int>     <dbl>
#> 1      1     0.549
#> 2      2     0.524
#> 3      3     0.477
```

### Efficient Block Handling

Vector-Based RSA automatically handles block structure:

``` r

# Compare with different block structures
block_sizes <- c(5, 10)
results_blocks <- lapply(block_sizes, function(size) {
  blocks <- rep(1:(n_stimuli/size), each = size)
  design <- vector_rsa_design(
    D = stim_distances,
    labels = stim_labels,
    block_var = blocks
  )
  model <- vector_rsa_model(
    dataset = dset,
    design = design,
    distfun = cordist()
  )
  suppressWarnings(suppressMessages(run_regional(model, region_mask)))
})

# Compare results
for (i in seq_along(block_sizes)) {
  cat("\nBlock size:", block_sizes[i], "\n")
  print(results_blocks[[i]]$performance_table)
}
#> 
#> Block size: 5 
#> # A tibble: 3 × 2
#>   roinum rsa_score
#>    <int>     <dbl>
#> 1      1     0.552
#> 2      2     0.522
#> 3      3     0.478
#> 
#> Block size: 10 
#> # A tibble: 3 × 2
#>   roinum rsa_score
#>    <int>     <dbl>
#> 1      1     0.554
#> 2      2     0.531
#> 3      3     0.482
```

## When to Use Each Method

### Feature‑Based RSA

Choose Feature‑Based RSA when you have a meaningful feature space and
care about predicting or reconstructing those dimensions from neural
activity (e.g., semantic features in language regions, low‑level visual
features in early visual cortex).

### Vector‑Based RSA

Use Vector‑Based RSA for large datasets or when block handling is
central to your design. Working with vectorized distances keeps memory
use in check and makes across‑block comparisons straightforward.

### Standard RSA

Standard RSA remains the best choice when your primary goal is to test
theoretical similarity structures and neither strict block handling nor
extreme memory efficiency is required.

## Summary

`rMVPA` offers three complementary RSA workflows: standard RSA for model
testing, Feature‑Based RSA for direct feature reconstruction, and
Vector‑Based RSA for efficient block‑aware distance comparisons. Pick
the approach that matches your research goal, data structure, and
computational budget.

For implementation details, see the reference pages for
[`feature_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_model.md),
[`vector_rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/vector_rsa_model.md),
and
[`rsa_model()`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md).
