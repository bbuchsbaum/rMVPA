# Validate an MVPA Analysis for Common Methodological Issues

Performs static checks on a model specification (or design + CV pair) to
detect potential data leakage, class imbalance, inadequate sample sizes,
and other common pitfalls *before* running a potentially expensive
analysis.

## Usage

``` r
validate_analysis(x, ...)

# S3 method for class 'mvpa_model'
validate_analysis(x, verbose = TRUE, ...)

# S3 method for class 'model_spec'
validate_analysis(x, verbose = TRUE, ...)

# S3 method for class 'mvpa_design'
validate_analysis(x, crossval = NULL, dataset = NULL, verbose = TRUE, ...)
```

## Arguments

- x:

  An object to validate. Typically an `mvpa_model` or an `mvpa_design`
  (in which case `crossval` must also be supplied).

- ...:

  Additional arguments (currently unused).

- verbose:

  Logical; if `TRUE` (the default), print the validation report to the
  console.

- crossval:

  A cross-validation specification (only needed when `x` is an
  `mvpa_design`).

- dataset:

  An optional `mvpa_dataset`, used for feature-count checks.

## Value

A `validation_result` object (invisibly) containing a list of individual
check results, each with `name`, `status` (`"pass"`, `"warn"`, or
`"fail"`), and `message`.

## Details

The following checks are performed:

- cv_block_alignment:

  Verifies that the CV scheme respects the blocking structure (e.g.
  run-level blocking for fMRI).

- block_structure:

  Checks that blocks are well-formed — enough blocks, roughly equal
  sizes, and all classes represented.

- class_balance:

  Checks class balance across folds and flags severe imbalance.

- single_class_folds:

  Detects folds where train or test sets contain only one class
  (guaranteed failure for classification).

- fold_sizes:

  Warns about very small test sets or extreme train/test ratios.

- sample_size:

  Warns when the total number of observations per class is very low.

- temporal_leakage:

  For non-blocked CV, detects when adjacent observations (likely
  temporally autocorrelated) are split across train and test.

## Examples

``` r
# Create a simple design
des_df <- data.frame(
  condition = factor(rep(c("A","B"), each = 50)),
  run = rep(1:5, each = 20)
)
des <- mvpa_design(des_df, y_train = ~ condition, block_var = ~ run)
cv  <- blocked_cross_validation(des$block_var)

# Validate — should pass all checks
validate_analysis(des, crossval = cv)
#> 
#>  MVPA Analysis Validation 
#> 
#>    [PASS]   cv_block_alignment 
#>     CV scheme is consistent with the blocking structure. 
#> 
#>    [PASS]   block_count 
#>     5 blocks found. 
#> 
#>    [WARN]   block_class_coverage 
#>     4 of 5 block(s) are missing at least one class. 
#> 
#>    [PASS]   min_observations 
#>     2 classes with 50 to 50 observations each. 
#> 
#>    [PASS]   class_balance 
#>     Class balance OK across folds (worst ratio: 1.7:1). 
#> 
#>    [WARN]   single_class_test 
#>     4 fold(s) have only one class in the test set (folds: 1, 2, 4,
#>     5). Metrics will be unreliable. 
#> 
#>    [PASS]   fold_sizes 
#>     Fold sizes OK (train: 80-80, test: 20-20). 
#> 
#>    Summary: 5 passed, 2 warnings 
#> 

# Risky: kfold CV ignoring run structure
cv_bad <- kfold_cross_validation(len = 100, nfolds = 5)
validate_analysis(des, crossval = cv_bad)
#> 
#>  MVPA Analysis Validation 
#> 
#>    [FAIL]   cv_block_alignment 
#>     Design has a blocking variable (e.g. fMRI run) but CV scheme is
#>     kfold_cross_validation, which ignores block structure.
#>     Observations from the same run may appear in both train and test,
#>     causing temporal autocorrelation leakage. Use
#>     blocked_cross_validation(block_var) instead. 
#> 
#>    [PASS]   block_count 
#>     5 blocks found. 
#> 
#>    [WARN]   block_class_coverage 
#>     4 of 5 block(s) are missing at least one class. 
#> 
#>    [PASS]   min_observations 
#>     2 classes with 50 to 50 observations each. 
#> 
#>    [PASS]   class_balance 
#>     Class balance OK across folds (worst ratio: 1.2:1). 
#> 
#>    [PASS]   single_class_folds 
#>     All folds have multiple classes in both train and test sets. 
#> 
#>    [PASS]   fold_sizes 
#>     Fold sizes OK (train: 80-80, test: 20-20). 
#> 
#>    Summary: 5 passed, 1 warnings, 1 failures 
#> 
```
