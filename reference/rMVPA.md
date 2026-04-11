# rMVPA: A package for multi-voxel pattern analysis (MVPA)

The rMVPA package provides a comprehensive suite of tools for advanced
neuroimaging data analysis using multi-voxel pattern analysis (MVPA). It
supports various techniques including region-of-interest (ROI) and
searchlight analyses. Key functionalities cover:

- Classification-based MVPA

- Representational Similarity Analysis (RSA), including standard RSA,
  vector-based RSA, and feature-based RSA.

- Contrast RSA (MS-ReVE style analyses)

- Flexible cross-validation schemes

- Feature selection methods

- Tools for constructing and managing MVPA datasets and designs.

## See also

Useful functions:

- [`mvpa_dataset`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
  for creating datasets

- [`run_searchlight`](http://bbuchsbaum.github.io/rMVPA/reference/run_searchlight.md)
  for searchlight analyses

- [`run_regional`](http://bbuchsbaum.github.io/rMVPA/reference/run_regional-methods.md)
  for ROI analyses

- [`spatial_nmf_maps`](http://bbuchsbaum.github.io/rMVPA/reference/spatial_nmf_maps.md)
  for spatial NMF and component inference

- [`rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md),
  [`contrast_rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/contrast_rsa_model.md),
  [`vector_rsa_model`](http://bbuchsbaum.github.io/rMVPA/reference/vector_rsa_model.md)
  for different RSA approaches

For a detailed introduction and examples, see the package vignette:
`vignette("rMVPA_introduction", package = "rMVPA")` (You may need to
build vignettes first if installing from source).

## Author

**Maintainer**: Bradley Buchsbaum <brad.buchsbaum@gmail.com>
([ORCID](https://orcid.org/0000-0002-1108-4866))
