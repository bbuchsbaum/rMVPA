# Model Targets

Retrieve the model-specific targets of a design together with their
meaning. Cross-validation folds are built from `cv_labels`; the
estimator sees the targets returned here. The two are identical for
ordinary classification designs but differ when a design carries
matrix-valued or continuous targets (feature prediction, feature sets,
feature RSA).

## Usage

``` r
model_targets(design, partition = c("train", "test"), ...)

# S3 method for class 'mvpa_design'
model_targets(design, partition = c("train", "test"), ...)

# S3 method for class 'feature_rsa_design'
model_targets(design, partition = c("train", "test"), ...)

# S3 method for class 'feature_sets_design'
model_targets(design, partition = c("train", "test"), ...)
```

## Arguments

- design:

  A design object such as
  [`mvpa_design`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md),
  [`feature_sets_design`](https://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md),
  or
  [`feature_rsa_design`](https://bbuchsbaum.github.io/rMVPA/reference/feature_rsa_design.md).

- partition:

  Either `"train"` or `"test"`. Requesting the test partition of a
  design without test targets returns `NULL`.

- ...:

  Additional arguments passed to methods.

## Value

A `model_targets` list as described above, or `NULL` when the requested
partition has no targets.

## Details

`model_targets` completes the `cv_labels` / `targets` split introduced
in
[`mvpa_design`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md).
It never changes what
[`y_train`](https://bbuchsbaum.github.io/rMVPA/reference/y_train-methods.md)
returns:
[`y_train()`](https://bbuchsbaum.github.io/rMVPA/reference/y_train-methods.md)
keeps returning the cross-validation labels for backward compatibility.

The returned object is a list of class `model_targets` with:

- values:

  A factor, numeric vector, or numeric matrix with one row (or element)
  per observation of the requested partition.

- observation_ids:

  Integer row identifiers, aligned with `values`.

- response_ids:

  Character identifiers for the responses: factor levels, matrix column
  names, or `NULL` for an unnamed scalar response.

- response_groups:

  Optional grouping of responses (for example the feature-set membership
  of a
  [`feature_sets_design`](https://bbuchsbaum.github.io/rMVPA/reference/feature_sets_design.md)),
  or `NULL`.

- row_weights:

  Optional numeric observation weights, or `NULL`.

- type:

  One of `"categorical"`, `"continuous"`, or `"matrix"`.

- partition:

  The requested partition.

## See also

[`mvpa_design`](https://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md),
[`y_train`](https://bbuchsbaum.github.io/rMVPA/reference/y_train-methods.md)

## Examples

``` r
des <- mvpa_design(data.frame(cond = rep(c("a", "b"), 10)), y_train = ~ cond)
model_targets(des)$type
#> [1] "categorical"

feats <- matrix(rnorm(20 * 3), 20, 3, dimnames = list(NULL, c("f1", "f2", "f3")))
des2 <- mvpa_design(data.frame(id = 1:20), cv_labels = 1:20, targets = feats)
dim(model_targets(des2)$values)
#> [1] 20  3
```
