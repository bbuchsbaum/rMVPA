# Create a Distance Function Object

Constructs a generic distance function object, storing:

- `name`: The name (method) of the distance (e.g., "euclidean",
  "cordist", "mahalanobis").

- `labels`: (Optional) a vector of labels associated with the data rows.

- `...`: Additional parameters relevant to the specific distance method
  (e.g., correlation method for `cordist`, number of components for
  `pcadist`, etc.).

## Usage

``` r
create_dist(name, labels = NULL, center = c("none", "stimulus_mean"), ...)
```

## Arguments

- name:

  A character string specifying the distance method (e.g., "euclidean",
  "cordist").

- labels:

  A vector of row labels (optional), primarily for
  informational/reference purposes.

- center:

  Optional pattern-centering method applied to `X` before distances are
  computed. Use `"stimulus_mean"` to subtract the across-stimulus mean
  pattern (Hanson-style).

- ...:

  Additional parameters for the distance method (e.g.
  \`method="pearson"\` for correlation, or `whiten=TRUE` for PCA-based
  distances).

## Value

An S3 object with class `c(name, "distfun")` that can be passed to
`pairwise_dist()`.

## Details

This object is used by `pairwise_dist()` to compute an **N x N** matrix
of pairwise distances between rows of a data matrix.

The distance function object itself does *not* exclude same-block
comparisons or reorder rows by label. Those tasks (if needed) are
handled downstream (for example, in `second_order_similarity`).

## Examples

``` r
# Create a Euclidean distance function object
dist_obj_euc <- create_dist("euclidean")

# Create a correlation distance function object with a specified correlation method
dist_obj_cor <- create_dist("cordist", method="spearman")
```
