# Balance Cross-Validation Partitions

Modifies a cross-validation partitioning scheme to ensure that each
target class has an equal number of samples within each training fold,
and optionally within each test fold, using either sub-sampling or
oversampling.

## Usage

``` r
balance_partitions(obj, design, ...)

# Default S3 method
balance_partitions(obj, design, method = "subsample", ...)

# S3 method for class 'blocked_cross_validation'
balance_partitions(obj, design, method = "subsample", ...)

# S3 method for class 'kfold_cross_validation'
balance_partitions(obj, design, method = "subsample", ...)

# S3 method for class 'twofold_blocked_cross_validation'
balance_partitions(obj, design, method = "subsample", ...)

# S3 method for class 'bootstrap_blocked_cross_validation'
balance_partitions(obj, design, method = "subsample", ...)

# S3 method for class 'sequential_blocked_cross_validation'
balance_partitions(obj, design, method = "subsample", ...)

# S3 method for class 'custom_cross_validation'
balance_partitions(obj, design, method = "subsample", ...)
```

## Arguments

- obj:

  A \`cross_validation\` object (e.g., from
  \`blocked_cross_validation\`, \`kfold_cross_validation\`).

- design:

  An \`mvpa_design\` object containing the target labels
  (\`.sa.targets\`) corresponding to the original dataset samples.

- ...:

  Additional arguments passed to internal balancing functions:

  balance_test_set

  :   Logical. If \`TRUE\` (default), balance the test sets in each fold
      as well using the specified \`method\`. \*\*Note:\*\* Oversampling
      the test set is generally not recommended as it can lead to
      misleading performance estimates. A warning will be issued if
      \`balance_test_set=TRUE\` and \`method="oversample"\`.

  seed

  :   An optional integer seed for the random number generator for
      reproducible sampling. If \`NULL\` (default), the result varies.

- method:

  Character string specifying the balancing method: - \`"subsample"\`
  (default): Down-samples majority classes to match the size of the
  smallest class (sampling without replacement). - \`"oversample"\`:
  Up-samples minority classes to match the size of the largest class
  (sampling with replacement).

## Value

A \`custom_cross_validation\` object where the sample indices in
\`.train_indices\` (and optionally \`.test_indices\`) for each fold have
been resampled to ensure balanced representation of target classes.

## Details

\*\*Sub-sampling (\`method="subsample"\`)\*\*: 1. Identifies the target
class with the minimum number of samples (\`min_count\`) within the set
(train or test). 2. For \*every\* target class within that set, it
randomly selects exactly \`min_count\` samples \*without replacement\*.
3. Discards samples from majority classes.

\*\*Oversampling (\`method="oversample"\`)\*\*: 1. Identifies the target
class with the maximum number of samples (\`max_count\`) within the set
(train or test). 2. For \*every\* target class within that set, it
randomly selects exactly \`max_count\` samples \*with replacement\*. 3.
Duplicates samples from minority classes.

Balancing guarantees that after the process, each target class appears
equally often within each balanced training (and optionally testing)
set. This is useful for preventing classifiers from being biased towards
majority classes.

The output is always a \`custom_cross_validation\` object because the
balancing process defines specific, explicit sets of indices for each
fold, which may no longer strictly adhere to the original blocking or
k-fold structure.

## Examples

``` r
# Create an imbalanced dataset design (more class 'b')
design_df <- data.frame(condition = factor(rep(c("a", "b", "b"), 20)),
                       block = rep(1:6, each = 10))
des <- mvpa_design(design_df, y_train = ~ condition, block_var = ~ block)

# Create standard blocked partitions (likely unbalanced)
cval_unbalanced <- blocked_cross_validation(des$block_var)
print("Unbalanced Counts (Example Fold 1 Train):")
#> [1] "Unbalanced Counts (Example Fold 1 Train):"
print(table(des$y_train[unlist(crossval_samples(cval_unbalanced,
          design_df, des$y_train)$train[[1]]$idx)]))
#> 
#>  a  b 
#> 16 34 

# Balance using sub-sampling (default)
cval_sub <- balance_partitions(cval_unbalanced, des, seed = 1)
print(cval_sub)
#> 
#> Custom Cross-Validation
#> 
#> - Configuration
#>   - Observations:    60
#>   - Number of Folds: 6
#> 
#> - Folds
#>   Fold    Train    Test   
#>   1       32       8      
#>   2       34       6      
#>   3       34       6      
#>   4       32       8      
#>   5       34       6      
#>   6       34       6      
#>   (mean train = 33, mean test = 7)
#> 
print("Subsample Balanced Counts (Example Fold 1 Train):")
#> [1] "Subsample Balanced Counts (Example Fold 1 Train):"
print(table(crossval_samples(cval_sub, design_df, des$y_train)$ytrain[[1]]))
#> 
#>  a  b 
#> 16 16 

# Balance using over-sampling
cval_over <- balance_partitions(cval_unbalanced, des, method = "oversample", seed = 2)
#> Warning: Oversampling the test set ('balance_test_set = TRUE', method = 'oversample') is generally not recommended and may lead to inflated performance metrics.
print(cval_over)
#> 
#> Custom Cross-Validation
#> 
#> - Configuration
#>   - Observations:    60
#>   - Number of Folds: 6
#> 
#> - Folds
#>   Fold    Train    Test   
#>   1       68       12     
#>   2       66       14     
#>   3       66       14     
#>   4       68       12     
#>   5       66       14     
#>   6       66       14     
#>   (mean train = 67, mean test = 13)
#> 
print("Oversample Balanced Counts (Example Fold 1 Train):")
#> [1] "Oversample Balanced Counts (Example Fold 1 Train):"
print(table(crossval_samples(cval_over, design_df, des$y_train)$ytrain[[1]]))
#> 
#>  a  b 
#> 34 34 
```
