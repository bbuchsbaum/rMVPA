# Distance Function Constructors

These convenience functions build specific types of distance function
objects via
[`create_dist`](http://bbuchsbaum.github.io/rMVPA/reference/create_dist.md).
Each returns an S3 object inheriting from `c("<method>", "distfun")`.

## Usage

``` r
cordist(
  labels = NULL,
  method = c("pearson", "spearman"),
  center = c("none", "stimulus_mean")
)

mahadist(labels = NULL, center = c("none", "stimulus_mean"))

eucdist(labels = NULL, center = c("none", "stimulus_mean"))

euclidean(labels = NULL, center = c("none", "stimulus_mean"))

robustmahadist(labels = NULL, center = c("none", "stimulus_mean"))

pcadist(
  labels = NULL,
  ncomp = 2,
  whiten = TRUE,
  threshfun = NULL,
  center = c("none", "stimulus_mean"),
  dist_method = c("euclidean", "manhattan", "cosine")
)
```

## Arguments

- labels:

  Optional vector of row labels (not directly used in distance
  calculation).

- method:

  For `cordist`, the correlation method: "pearson" or "spearman".

- center:

  Optional pattern-centering method applied to `X` before distances are
  computed. Use `"stimulus_mean"` to subtract the across-stimulus mean
  pattern (Hanson-style).

- ncomp:

  For `pcadist`, the number of components (or a function threshold).

- whiten:

  For `pcadist`, whether to whiten principal components (logical).

- threshfun:

  For `pcadist`, an optional function that determines how many PCs to
  retain based on `pres$sdev^2`.

- dist_method:

  For `pcadist`, the base distance measure in PC space ("euclidean",
  "manhattan", or "cosine").

## Value

An S3 object with class `c("<method>", "distfun")`.

## Details

\- `cordist(labels, method="pearson")` -\> correlation-based distance. -
`mahadist(labels)` -\> Mahalanobis distance. - `eucdist(labels)` -\>
Euclidean distance. - `robustmahadist(labels)` -\> Mahalanobis distance
using robust covariance. - `pcadist(labels, ...)` -\> distance in
reduced PCA space.

## See also

[`create_dist`](http://bbuchsbaum.github.io/rMVPA/reference/create_dist.md)
for the underlying constructor.

## Examples

``` r
dist_obj_1 <- cordist(method="pearson")
dist_obj_2 <- mahadist()
dist_obj_3 <- eucdist()
dist_obj_4 <- robustmahadist()
dist_obj_5 <- pcadist(ncomp=2, dist_method="cosine")
```
