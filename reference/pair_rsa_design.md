# Construct a pair-observation RSA design

Generalizes
[`rsa_design`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md)
to support arbitrary pair-observation geometries: lower-triangle
within-domain pairs (the classical RSA layout), rectangular
between-domain pairs (e.g. items in domain A vs. items in domain B), and
function-valued model RDM entries. Function entries may either operate
on item identifiers, or on item identifiers plus feature rows when
`features_a` / `features_b` are supplied.

## Usage

``` r
pair_rsa_design(
  items_a,
  items_b = NULL,
  model = list(),
  nuisance = list(),
  pairs = c("within", "between"),
  features_a = NULL,
  features_b = NULL,
  block_var_a = NULL,
  block_var_b = NULL,
  keep_intra_run = FALSE,
  row_idx_a = NULL,
  row_idx_b = NULL
)
```

## Arguments

- items_a:

  Vector of item identifiers in domain A. Length defines `n_a`.

- items_b:

  Optional vector of item identifiers in domain B. Required when
  `pairs = "between"`; ignored otherwise. Length defines `n_b`.

- model:

  Named list of model RDM specifications. Each entry may be a `dist`
  object, a numeric square matrix (within mode) or `n_a x n_b` matrix
  (between mode), a numeric vector of length `n_pairs`, or a function
  `function(a, b)` returning a numeric vector of pairwise values for the
  requested pair table. If `features_a` is supplied and the function has
  at least four formal arguments (or `...`), it is called as
  `function(a, b, features_a, features_b)` where the feature arguments
  are row-aligned pair tables.

- nuisance:

  Optional named list of nuisance pair predictors using the same
  accepted forms as `model`. Included in the RSA design matrix but
  excluded from model-space fingerprints returned by
  `rsa_model(..., return_fingerprint = TRUE)`.

- pairs:

  Either `"within"` (default) or `"between"`.

- features_a:

  Optional data frame or matrix of item features for `items_a`.
  Function-valued model/nuisance entries can use these rows to define
  feature-pair dissimilarities.

- features_b:

  Optional data frame or matrix of item features for `items_b`. Defaults
  to `features_a` in within mode.

- block_var_a:

  Optional vector of block labels (length `n_a`). If supplied and
  `keep_intra_run = FALSE`, within-block pairs are excluded.

- block_var_b:

  Optional vector of block labels for domain B, length `n_b`. Defaults
  to `block_var_a` in within mode.

- keep_intra_run:

  Logical; if `TRUE`, do not drop within-block pairs.

- row_idx_a:

  Optional integer vector of dataset row indices corresponding to
  `items_a`. In within mode, this lets a pair design address a
  subset/reordering of dataset rows. Required for `pairs = "between"` so
  the per-ROI engine can extract the correct neural sub-blocks.

- row_idx_b:

  Integer vector of dataset row indices for `items_b`. Required for
  `pairs = "between"`.

## Value

A list with class `c("pair_rsa_design", "rsa_design", "list")`
containing all fields produced by
[`rsa_design()`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md)
plus

- pair_kind:

  Either `"within"` or `"between"`.

- items_a, items_b:

  Item identifier vectors.

- n_a, n_b:

  Item counts.

- pair_index:

  A data frame describing each retained pair.

- row_idx_a, row_idx_b:

  Dataset row indices when `pairs = "between"`.

## Details

The returned object inherits from `rsa_design` so it is a drop-in
replacement when used with
[`rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md)
and the regional / searchlight engines. Within-domain pair designs are
fully interoperable with the existing `train_model.rsa_model` path;
between-domain designs are dispatched on the `pair_kind` field, which
causes `train_model.rsa_model` to compute a rectangular neural-pair
dissimilarity block instead of the lower triangle.

## See also

[`rsa_design`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_design.md),
[`rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md),
[`model_space_connectivity`](http://bbuchsbaum.github.io/rMVPA/reference/model_space_connectivity.md)

## Examples

``` r
set.seed(1)
items <- paste0("item", 1:8)
R1 <- as.matrix(dist(matrix(rnorm(8 * 4), 8, 4)))
R2 <- as.matrix(dist(matrix(rnorm(8 * 4), 8, 4)))
rownames(R1) <- colnames(R1) <- rownames(R2) <- colnames(R2) <- items
des <- pair_rsa_design(items, model = list(rdm1 = R1, rdm2 = R2))
lengths(des$model_mat)
#> rdm1 rdm2 
#>   28   28 
```
