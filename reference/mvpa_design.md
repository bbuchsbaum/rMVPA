# Create an MVPA Design Object

Creates a design object for MVPA analysis that encapsulates training and
testing designs, response variables, and optional blocking and splitting
factors.

## Usage

``` r
mvpa_design(
  train_design,
  test_design = NULL,
  y_train = NULL,
  y_test = NULL,
  block_var = NULL,
  split_by = NULL,
  cv_labels = NULL,
  targets = NULL,
  ...
)
```

## Arguments

- train_design:

  A data frame containing the training design matrix

- test_design:

  Optional data frame containing the test design matrix (default: NULL)

- y_train:

  Formula or vector specifying the training response variable (old path;
  mutually exclusive with `cv_labels`)

- y_test:

  Optional formula or vector specifying the test response variable
  (default: NULL)

- block_var:

  Optional formula or vector specifying the blocking variable for
  cross-validation

- split_by:

  Optional formula or vector for splitting analyses

- cv_labels:

  Optional vector of labels used for cross-validation fold construction
  (new path; mutually exclusive with `y_train`). If provided without
  `targets`, `targets` defaults to `cv_labels`.

- targets:

  Optional vector of model-specific training targets. Defaults to
  `cv_labels` when `cv_labels` is supplied, or to the parsed `y_train`
  when the old path is used.

- ...:

  Additional arguments (currently unused)

## Value

An `mvpa_design` object (S3 class) containing:

- train_design:

  Data frame of training design

- test_design:

  Data frame of test design (if provided)

- cv_labels:

  Labels used for cross-validation fold construction

- targets:

  Model-specific training targets

- y_train:

  Alias for `cv_labels` (for backward compatibility)

- y_test:

  Test response variable (if provided)

- block_var:

  Blocking variable for cross-validation (if provided)

- split_by:

  Splitting factor (if provided)

## Details

The `y_train` and `y_test` can be specified either as formulas (e.g., ~
condition) or as vectors. If formulas are used, they are evaluated
within the respective design matrices.

The `block_var` and `split_by` can also be specified as formulas or
vectors. If formulas, they are evaluated within the training design
matrix.

The new `cv_labels` and `targets` parameters allow separating the labels
used for fold construction from the actual training targets. When using
the old `y_train` path, both `cv_labels` and `targets` are set to the
parsed `y_train` value.

## See also

[`mvpa_dataset`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
for creating the corresponding dataset object

## Examples

``` r
# Basic design with only training data
train_df <- data.frame(condition = rep(c("A", "B"), each = 50),
                       block = rep(1:5, each = 20),
                       group = rep(c("Group1", "Group2"), 50))
design <- mvpa_design(train_df, y_train = ~ condition)

# Design with test data and blocking variable
test_df <- data.frame(condition = rep(c("A", "B"), each = 25))
design_with_test <- mvpa_design(
  train_df,
  test_df,
  y_train = ~ condition,
  y_test = ~ condition,
  block_var = ~ block
)

# Design with split_by factor
design_split <- mvpa_design(
  train_df,
  y_train = ~ condition,
  split_by = ~ group
)
```
