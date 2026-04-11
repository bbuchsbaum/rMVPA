# Cross-Validation Strategies in rMVPA

## Introduction

Cross-validation is a critical component in evaluating the performance
and generalizability of MVPA models. The `rMVPA` package provides
several cross-validation strategies specifically designed for
neuroimaging data, where temporal structure and run/block organization
must be carefully considered.

This vignette explains why cross‑validation matters for MVPA, introduces
schemes tailored to neuroimaging data, shows how to implement them, and
closes with practical guidance and examples.

## Cross-Validation Fundamentals

### Why Cross-Validation Matters in MVPA

Cross-validation is essential in neuroimaging analyses. It gives us
unbiased estimates of how well our models perform and shows whether
brain activity patterns truly generalize to new data. By separating
training and test sets, cross-validation prevents us from overfitting to
noise in the data. It also maintains temporal independence between sets,
which is crucial for time-series neuroimaging data.

### Available cross‑validation schemes

`rMVPA` supports blocked CV (using runs as folds), k‑fold CV (random
partitions), bootstrap blocked CV (resampling within runs), sequential
blocked CV (ordered folds within runs), two‑fold blocked CV (fast split
across runs), and fully custom schemes when you need explicit train/test
indices.

## Blocked Cross-Validation

### Overview

Blocked cross-validation is the most common approach for fMRI data. It
respects the temporal structure of the data by using scanning runs as
natural validation blocks.

``` r
# Create a simple blocked structure: 5 runs with 20 trials each
block_var <- rep(1:5, each = 20)
cval <- blocked_cross_validation(block_var)
print(cval)
```

    ## 
    ##  Blocked Cross-Validation 
    ## 
    ## - Dataset Information 
    ##   - Observations:  100 
    ##   - Number of Folds:  5 
    ## - Block Information 
    ##   - Total Blocks:  5 
    ##   - Mean Block Size:  20  (SD:  0 ) 
    ##   - Block Sizes:  1: 20, 2: 20, 3: 20, 4: 20, 5: 20

### Implementation Example

``` r
# Generate example data
set.seed(123)
dat <- data.frame(
  x1 = rnorm(100),  # 100 trials total
  x2 = rnorm(100),
  x3 = rnorm(100)
)
y <- factor(rep(letters[1:5], length.out = 100))  # 5 conditions

# Generate cross-validation samples
samples <- crossval_samples(cval, dat, y)
print(samples)
```

    ## # A tibble: 5 × 5
    ##   ytrain       ytest        train               test                .id  
    ##   <named list> <named list> <named list>        <named list>        <chr>
    ## 1 <fct [80]>   <fct [20]>   <resample [80 x 3]> <resample [20 x 3]> 01   
    ## 2 <fct [80]>   <fct [20]>   <resample [80 x 3]> <resample [20 x 3]> 02   
    ## 3 <fct [80]>   <fct [20]>   <resample [80 x 3]> <resample [20 x 3]> 03   
    ## 4 <fct [80]>   <fct [20]>   <resample [80 x 3]> <resample [20 x 3]> 04   
    ## 5 <fct [80]>   <fct [20]>   <resample [80 x 3]> <resample [20 x 3]> 05

### Understanding the output

Each row in the samples tibble contains the training and test labels
(`ytrain`, `ytest`), the corresponding data subsets (`train`, `test`),
and a fold identifier (`.id`).

## Bootstrap Blocked Cross-Validation

### Overview

This method combines blocking with bootstrap resampling, providing more
stable performance estimates while respecting the run structure.

``` r
# Create bootstrap blocked CV with 20 repetitions
boot_cval <- bootstrap_blocked_cross_validation(block_var, nreps = 20)
print(boot_cval)
```

    ## 
    ##  Bootstrap Blocked Cross-Validation 
    ## 
    ## - Configuration 
    ##   - Observations:  100 
    ##   - Bootstrap Repetitions:  20 
    ## - Block Information 
    ##   - Total Blocks:  5 
    ##   - Mean Block Size:  20  (SD:  0 ) 
    ##   - Block Sizes:  1: 20, 2: 20, 3: 20, 4: 20, 5: 20 
    ## - Sampling Weights 
    ##   - Status:  None  (uniform sampling)

``` r
# Generate samples
boot_samples <- crossval_samples(boot_cval, dat, y)
print(boot_samples)
```

    ## # A tibble: 100 × 5
    ##    ytrain     ytest      train               test                .id  
    ##    <list>     <list>     <list>              <list>              <chr>
    ##  1 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 001  
    ##  2 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 002  
    ##  3 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 003  
    ##  4 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 004  
    ##  5 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 005  
    ##  6 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 006  
    ##  7 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 007  
    ##  8 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 008  
    ##  9 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 009  
    ## 10 <fct [80]> <fct [20]> <resample [80 x 3]> <resample [20 x 3]> 010  
    ## # ℹ 90 more rows

### Optional Weighted Sampling

You can provide weights to influence the sampling probability:

``` r
# Create weights (e.g., based on motion parameters)
weights <- runif(length(block_var))
weighted_boot_cval <- bootstrap_blocked_cross_validation(
  block_var, 
  nreps = 20,
  weights = weights
)
print(weighted_boot_cval)
```

    ## 
    ##  Bootstrap Blocked Cross-Validation 
    ## 
    ## - Configuration 
    ##   - Observations:  100 
    ##   - Bootstrap Repetitions:  20 
    ## - Block Information 
    ##   - Total Blocks:  5 
    ##   - Mean Block Size:  20  (SD:  0 ) 
    ##   - Block Sizes:  1: 20, 2: 20, 3: 20, 4: 20, 5: 20 
    ## - Sampling Weights 
    ##   - Status:  Present 
    ##   - Range:  [0.001, 0.020] 
    ##   - Non-zero Weights:  100  ( 100.0% )

## Sequential Blocked Cross-Validation

### Overview

This method creates sequential folds within each block, useful when
temporal order matters.

``` r
# Create sequential blocked CV with 2 folds and 4 repetitions
seq_cval <- sequential_blocked_cross_validation(
  block_var,
  nfolds = 2,
  nreps = 4
)
print(seq_cval)
```

    ## $block_var
    ##   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
    ##  [38] 2 2 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 4 4 4 4 4 4 4 4 4 4 4 4 4 4
    ##  [75] 4 4 4 4 4 4 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5
    ## 
    ## $nfolds
    ## [1] 2
    ## 
    ## $nreps
    ## [1] 4
    ## 
    ## $block_ind
    ## [1] 1 2 3 4 5
    ## 
    ## attr(,"class")
    ## [1] "sequential_blocked_cross_validation" "cross_validation"                   
    ## [3] "list"

## Two-Fold Blocked Cross-Validation

### Overview

This approach randomly splits blocks into two groups, useful for rapid
performance estimation.

``` r
# Create two-fold blocked CV with 10 repetitions
twofold_cval <- twofold_blocked_cross_validation(block_var, nreps = 10)
print(twofold_cval)
```

    ## 
    ##  Two-Fold Blocked Cross-Validation 
    ## 
    ## - Configuration 
    ##   - Observations:  100 
    ##   - Number of Folds:  2 
    ##   - Repetitions:  10 
    ## - Block Information 
    ##   - Total Blocks:  5 
    ##   - Mean Block Size:  20  (SD:  0 ) 
    ##   - Block Sizes:  1: 20, 2: 20, 3: 20, 4: 20, 5: 20

## Custom Cross-Validation

### Overview

For specialized validation schemes, you can define custom
training/testing splits:

``` r
# Define custom splits
custom_splits <- list(
  list(train = 1:60, test = 61:100),
  list(train = 1:40, test = 41:100),
  list(train = 1:80, test = 81:100)
)

# Create custom CV
custom_cval <- custom_cross_validation(custom_splits)
print(custom_cval)
```

    ## $sample_set
    ## $sample_set[[1]]
    ## $sample_set[[1]]$train
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
    ## [51] 51 52 53 54 55 56 57 58 59 60
    ## 
    ## $sample_set[[1]]$test
    ##  [1]  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78  79
    ## [20]  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98
    ## [39]  99 100
    ## 
    ## 
    ## $sample_set[[2]]
    ## $sample_set[[2]]$train
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
    ## 
    ## $sample_set[[2]]$test
    ##  [1]  41  42  43  44  45  46  47  48  49  50  51  52  53  54  55  56  57  58  59
    ## [20]  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77  78
    ## [39]  79  80  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97
    ## [58]  98  99 100
    ## 
    ## 
    ## $sample_set[[3]]
    ## $sample_set[[3]]$train
    ##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
    ## [26] 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50
    ## [51] 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75
    ## [76] 76 77 78 79 80
    ## 
    ## $sample_set[[3]]$test
    ##  [1]  81  82  83  84  85  86  87  88  89  90  91  92  93  94  95  96  97  98  99
    ## [20] 100
    ## 
    ## 
    ## 
    ## $nfolds
    ## [1] 3
    ## 
    ## attr(,"class")
    ## [1] "custom_cross_validation" "cross_validation"       
    ## [3] "list"

## Practical Example: Model Training

Here’s a complete example using blocked cross-validation with an SDA
classifier:

``` r
# Setup cross-validation
block_var <- rep(1:5, each = 20)
cval <- blocked_cross_validation(block_var)

# Generate data
set.seed(123)
dat <- data.frame(matrix(rnorm(100 * 10), 100, 10))
y <- factor(rep(letters[1:5], 20))

# Generate CV samples
samples <- crossval_samples(cval, dat, y)

# Train models for each fold
model_fits <- samples %>%
  rowwise() %>%
  mutate(fit = list(sda::sda(as.matrix(as.data.frame(train)), ytrain, verbose = FALSE)))

print(model_fits)
```

    ## # A tibble: 5 × 6
    ## # Rowwise: 
    ##   ytrain       ytest       train                test                 .id   fit  
    ##   <named list> <named lis> <named list>         <named list>         <chr> <lis>
    ## 1 <fct [80]>   <fct [20]>  <resample [80 x 10]> <resample [20 x 10]> 01    <sda>
    ## 2 <fct [80]>   <fct [20]>  <resample [80 x 10]> <resample [20 x 10]> 02    <sda>
    ## 3 <fct [80]>   <fct [20]>  <resample [80 x 10]> <resample [20 x 10]> 03    <sda>
    ## 4 <fct [80]>   <fct [20]>  <resample [80 x 10]> <resample [20 x 10]> 04    <sda>
    ## 5 <fct [80]>   <fct [20]>  <resample [80 x 10]> <resample [20 x 10]> 05    <sda>

## Best practices

When choosing a cross‑validation strategy, align the scheme with your
data structure and analysis goals. Blocked CV should be your default for
fMRI because it respects run boundaries and temporal dependence.
Sequential CV is helpful when order matters. Bootstrap variants trade
additional computation for more stable estimates. Two‑fold CV offers a
quick sanity check when you need speed.

In small datasets, bootstrap blocked CV can stabilize estimates; in
larger datasets, simple blocked CV often suffices. Always keep train and
test data from the same run separate to avoid temporal leakage. Balance
stability against computation time: bootstrap and sequential methods are
costlier, whereas two‑fold CV provides fast approximations for early
exploration.

## Summary

Choose the right cross-validation strategy for your neuroimaging data.
Match the CV approach to your data structure, account for temporal
dependencies in fMRI, and use bootstrap methods when you need more
stable estimates. The rMVPA package supports all these approaches,
including custom schemes for specialized needs.

For implementation details, see
[`blocked_cross_validation()`](http://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md),
[`bootstrap_blocked_cross_validation()`](http://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md),
and
[`custom_cross_validation()`](http://bbuchsbaum.github.io/rMVPA/reference/cross_validation.md).

## Integration with Regional and Searchlight Analyses

Cross-validation strategies in `rMVPA` are designed to work seamlessly
with both regional and searchlight analyses. Here we’ll demonstrate how
to incorporate different cross-validation schemes into these analyses.

### Regional Analysis Example

Let’s perform a regional MVPA analysis using different cross-validation
strategies:

``` r
# Generate a sample dataset
data_out <- gen_sample_dataset(D = c(6,6,6), nobs = 80, blocks = 4, nlevels = 2)

# Create a region mask
mask <- data_out$dataset$mask
nvox <- sum(mask)
region_mask <- neuroim2::NeuroVol(
  sample(1:3, size = nvox, replace = TRUE), 
  space(mask), 
  indices = which(mask > 0)
)

# Create MVPA dataset
dset <- mvpa_dataset(data_out$dataset$train_data, mask = data_out$dataset$mask)

# Load the classification model
mod <- load_model("sda_notune")
tune_grid <- data.frame(lambda = 0.01, diagonal = FALSE)

# Example 1: Using Blocked Cross-Validation
blocked_cv <- blocked_cross_validation(data_out$design$block_var)
mvpa_mod_blocked <- mvpa_model(
  mod, 
  dataset = dset, 
  design = data_out$design,
  crossval = blocked_cv,
  tune_grid = tune_grid
)

# Run regional analysis with blocked CV
results_blocked <- run_regional(mvpa_mod_blocked, region_mask)

# Example 2: Using Bootstrap Blocked Cross-Validation
bootstrap_cv <- bootstrap_blocked_cross_validation(
  data_out$design$block_var,
  nreps = 10
)
mvpa_mod_boot <- mvpa_model(
  mod,
  dataset = dset,
  design = data_out$design,
  crossval = bootstrap_cv,
  tune_grid = tune_grid
)

# Run regional analysis with bootstrap CV
results_boot <- run_regional(mvpa_mod_boot, region_mask)

# Compare performance between CV strategies
cat("Blocked CV Performance:\n")
```

    ## Blocked CV Performance:

``` r
print(results_blocked$performance_table)
```

    ## # A tibble: 3 × 3
    ##   roinum Accuracy    AUC
    ##    <int>    <dbl>  <dbl>
    ## 1      1    0.562  0.191
    ## 2      2    0.6    0.305
    ## 3      3    0.462 -0.135

``` r
cat("\nBootstrap CV Performance:\n")
```

    ## 
    ## Bootstrap CV Performance:

``` r
print(results_boot$performance_table)
```

    ## # A tibble: 3 × 3
    ##   roinum Accuracy     AUC
    ##    <int>    <dbl>   <dbl>
    ## 1      1    0.562  0.176 
    ## 2      2    0.575  0.281 
    ## 3      3    0.475 -0.0988

### Searchlight Analysis Example

We can also use different cross-validation strategies in searchlight
analysis:

``` r
# Example 3: Searchlight with Sequential Blocked Cross-Validation
seq_cv <- sequential_blocked_cross_validation(
  data_out$design$block_var,
  nfolds = 2,
  nreps = 4
)

mvpa_mod_seq <- mvpa_model(
  mod,
  dataset = dset,
  design = data_out$design,
  crossval = seq_cv,
  tune_grid = tune_grid
)

# Run searchlight analysis
results_searchlight <- run_searchlight(
  mvpa_mod_seq,
  radius = 2,
  method = "standard"
)
```

### Key considerations

Searchlights multiply computation by the number of centers, so prefer
lighter schemes (two‑fold or blocked) for initial runs and reduce
`nreps` if using bootstrap. Use parallel processing where possible.
Expect minor differences in performance estimates across schemes;
bootstrap tends to be more conservative. For repeated schemes, average
over repetitions to summarize performance.

### Example: Comparing CV Strategies in Regional Analysis

Here’s a more detailed comparison of different CV strategies:

``` r
# Create different CV schemes
cv_schemes <- list(
  blocked = blocked_cross_validation(data_out$design$block_var),
  bootstrap = bootstrap_blocked_cross_validation(data_out$design$block_var, nreps = 10),
  twofold = twofold_blocked_cross_validation(data_out$design$block_var, nreps = 5)
)

# Run regional analysis with each CV scheme
results <- lapply(cv_schemes, function(cv) {
  mvpa_mod <- mvpa_model(
    mod,
    dataset = dset,
    design = data_out$design,
    crossval = cv,
    tune_grid = tune_grid
  )
  run_regional(mvpa_mod, region_mask)
})

# Compare performance across CV schemes
performance_comparison <- lapply(names(results), function(name) {
  perf <- results[[name]]$performance_table
  perf$cv_scheme <- name
  perf
})

# Combine results
all_performance <- do.call(rbind, performance_comparison)
print(all_performance)
```

    ## # A tibble: 9 × 4
    ##   roinum Accuracy     AUC cv_scheme
    ##    <int>    <dbl>   <dbl> <chr>    
    ## 1      1    0.562  0.191  blocked  
    ## 2      2    0.6    0.305  blocked  
    ## 3      3    0.462 -0.135  blocked  
    ## 4      1    0.562  0.231  bootstrap
    ## 5      2    0.55   0.206  bootstrap
    ## 6      3    0.538 -0.0400 bootstrap
    ## 7      1    0.562  0.114  twofold  
    ## 8      2    0.562  0.288  twofold  
    ## 9      3    0.462 -0.15   twofold

This integration demonstrates how different cross-validation strategies
can be easily incorporated into the broader MVPA analysis framework,
allowing researchers to choose the most appropriate validation approach
for their specific analysis needs.

For more details, see
[`vignette("Regional_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Regional_Analysis.md)
and
[`vignette("Searchlight_Analysis")`](http://bbuchsbaum.github.io/rMVPA/articles/Searchlight_Analysis.md).
