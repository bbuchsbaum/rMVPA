# Subsample Searchlight Centers

Selects a representative subset of searchlight centers for permutation
runs, optionally stratifying by the number of features per center.

## Usage

``` r
subsample_centers(
  dataset,
  searchlight,
  n_centers = NULL,
  fraction = 0.1,
  stratify_by = "nfeatures",
  redundancy_map = NULL,
  seed = NULL
)
```

## Arguments

- dataset:

  An `mvpa_dataset`.

- searchlight:

  A searchlight iterator (list of integer voxel-index vectors) as
  returned by
  [`get_searchlight()`](http://bbuchsbaum.github.io/rMVPA/reference/get_searchlight.md).

- n_centers:

  Optional integer count of centers to select. If `NULL`, derived from
  `fraction`.

- fraction:

  Numeric in (0, 1\]. Fraction of all centers to select when `n_centers`
  is `NULL`.

- stratify_by:

  Character. Covariate used for stratification (`"nfeatures"` only for
  now).

- redundancy_map:

  Optional named numeric vector of per-center redundancy values (not
  used unless `stratify_by = "redundancy"`).

- seed:

  Optional integer seed.

## Value

A list with:

- center_ids:

  Integer vector of selected center IDs.

- vox_list:

  Subsetted searchlight list.

- covariates:

  A `data.frame` with at least column `nfeatures`.
