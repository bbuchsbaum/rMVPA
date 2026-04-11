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
#> $sample_set
#> $sample_set[[1]]
#> $sample_set[[1]]$train
#>  [1] 11 12 13 14 16 17 19 22 23 24 25 27 28 30 31 32 34 37 39 40 41 43 46 47 49
#> [26] 50 52 55 56 58 59 60
#> 
#> $sample_set[[1]]$test
#> [1]  1  3  4  6  7  8  9 10
#> 
#> 
#> $sample_set[[2]]
#> $sample_set[[2]]$train
#>  [1]  1  2  4  6  7  9 10 21 22 23 25 27 28 29 31 32 34 36 37 38 39 40 43 44 45
#> [26] 46 48 49 50 52 53 55 56 58
#> 
#> $sample_set[[2]]$test
#> [1] 11 13 15 16 18 19
#> 
#> 
#> $sample_set[[3]]
#> $sample_set[[3]]$train
#>  [1]  1  2  3  4  7 10 11 12 13 15 16 17 18 19 31 33 34 35 37 39 40 42 43 46 49
#> [26] 50 51 52 53 55 57 58 59 60
#> 
#> $sample_set[[3]]$test
#> [1] 21 22 24 25 27 28
#> 
#> 
#> $sample_set[[4]]
#> $sample_set[[4]]$train
#>  [1]  1  2  3  4  6  7 10 11 13 16 18 19 22 23 24 25 28 29 30 41 43 45 46 47 48
#> [26] 49 50 51 52 55 58 59
#> 
#> $sample_set[[4]]$test
#> [1] 31 32 33 34 35 37 39 40
#> 
#> 
#> $sample_set[[5]]
#> $sample_set[[5]]$train
#>  [1]  1  2  3  4  5  7  8  9 10 13 14 16 18 19 21 22 23 24 25 26 28 30 31 34 37
#> [26] 40 51 52 53 54 55 57 58 59
#> 
#> $sample_set[[5]]$test
#> [1] 43 45 46 47 49 50
#> 
#> 
#> $sample_set[[6]]
#> $sample_set[[6]]$train
#>  [1]  1  3  4  6  7  8 10 11 12 13 14 15 16 19 21 22 25 26 28 30 31 32 34 37 39
#> [26] 40 41 42 43 44 45 46 48 49
#> 
#> $sample_set[[6]]$test
#> [1] 51 52 55 56 58 59
#> 
#> 
#> 
#> $nfolds
#> [1] 6
#> 
#> attr(,"class")
#> [1] "custom_cross_validation" "cross_validation"       
#> [3] "list"                   
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
#> $sample_set
#> $sample_set[[1]]
#> $sample_set[[1]]$train
#>  [1] 11 12 13 13 13 13 14 14 15 16 16 18 18 18 19 19 20 21 21 21 22 23 23 25 28
#> [26] 28 28 28 28 29 29 29 31 31 32 32 33 33 34 34 36 36 37 37 37 38 39 41 43 43
#> [51] 44 46 49 49 51 52 54 54 55 55 57 57 58 58 58 58 59 60
#> 
#> $sample_set[[1]]$test
#>  [1] 1 1 1 1 2 2 4 4 6 8 8 9
#> 
#> 
#> $sample_set[[2]]
#> $sample_set[[2]]$train
#>  [1]  1  1  1  3  4  5  6  7  7  7  8 10 22 22 22 22 23 23 23 25 25 28 28 29 29
#> [26] 30 32 32 34 34 34 35 36 36 37 38 38 38 40 40 40 41 41 42 43 43 44 45 46 47
#> [51] 49 49 49 50 51 52 52 54 54 55 56 56 56 58 59 60
#> 
#> $sample_set[[2]]$test
#>  [1] 11 13 14 16 16 16 18 18 18 19 19 19 20 20
#> 
#> 
#> $sample_set[[3]]
#> $sample_set[[3]]$train
#>  [1]  1  2  3  5  5  6  7  7  7  7  9 10 10 11 11 15 16 16 18 19 20 31 31 31 32
#> [26] 32 32 34 34 35 37 37 37 39 39 40 40 40 40 41 41 42 42 42 43 43 43 44 45 45
#> [51] 46 46 46 47 48 48 49 50 50 51 52 54 55 58 58 59
#> 
#> $sample_set[[3]]$test
#>  [1] 21 22 22 22 24 24 24 25 25 27 28 28 29 30
#> 
#> 
#> $sample_set[[4]]
#> $sample_set[[4]]$train
#>  [1]  1  3  6  7  7  7  9 10 10 12 13 13 13 14 14 14 14 16 17 17 18 19 19 21 21
#> [26] 23 24 25 25 25 25 27 27 27 28 28 28 42 42 43 43 43 45 45 45 45 46 46 46 46
#> [51] 46 48 48 49 50 52 52 53 53 53 54 55 55 55 57 57 58 59
#> 
#> $sample_set[[4]]$test
#>  [1] 31 32 33 33 34 35 37 37 37 38 39 40
#> 
#> 
#> $sample_set[[5]]
#> $sample_set[[5]]$train
#>  [1]  1  1  1  1  1  3  4  4  4  7  8  9  9 10 11 12 12 12 13 13 13 15 17 20 20
#> [26] 21 21 21 22 22 22 22 23 23 23 24 24 25 25 31 31 31 31 32 32 32 33 33 34 34
#> [51] 35 37 37 37 38 38 39 39 40 40 51 52 52 57 58 59
#> 
#> $sample_set[[5]]$test
#>  [1] 41 41 43 45 46 46 46 46 47 48 49 49 50 50
#> 
#> 
#> $sample_set[[6]]
#> $sample_set[[6]]$train
#>  [1]  1  3  5  7  8  8  9 11 11 13 13 13 14 14 15 15 16 16 16 17 18 18 19 20 21
#> [26] 21 22 24 25 25 25 26 26 27 28 28 28 28 30 30 31 31 32 33 34 34 34 39 39 40
#> [51] 40 41 42 42 43 43 44 46 46 47 48 49 49 49 49 49
#> 
#> $sample_set[[6]]$test
#>  [1] 52 52 53 54 54 54 55 55 57 57 58 58 58 59
#> 
#> 
#> 
#> $nfolds
#> [1] 6
#> 
#> attr(,"class")
#> [1] "custom_cross_validation" "cross_validation"       
#> [3] "list"                   
print("Oversample Balanced Counts (Example Fold 1 Train):")
#> [1] "Oversample Balanced Counts (Example Fold 1 Train):"
print(table(crossval_samples(cval_over, design_df, des$y_train)$ytrain[[1]]))
#> 
#>  a  b 
#> 34 34 
```
