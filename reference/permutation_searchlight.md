# Permutation-Based Inference for Searchlight MVPA

This module provides permutation testing infrastructure for searchlight
MVPA, including covariate-conditioned null distributions, stratified
subsampling, and FDR correction. It wraps existing rMVPA infrastructure;
a model type takes part by providing a
[`permute_labels`](https://bbuchsbaum.github.io/rMVPA/reference/permute_labels.md)
method for its design class.
