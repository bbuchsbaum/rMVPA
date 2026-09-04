# Permute Labels in a Design for Permutation Testing

Returns a copy of `design` in which the link between the brain data and
the quantity being decoded or modelled has been broken by a random
relabelling that respects the design's block structure. This is the
null-generating step of
[`run_permutation_searchlight`](https://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md).

## Usage

``` r
permute_labels(
  design,
  method = c("within_block", "circular_shift", "global"),
  seed = NULL
)

# Default S3 method
permute_labels(
  design,
  method = c("within_block", "circular_shift", "global"),
  seed = NULL
)

# S3 method for class 'mvpa_design'
permute_labels(
  design,
  method = c("within_block", "circular_shift", "global"),
  seed = NULL
)

# S3 method for class 'rsa_design'
permute_labels(
  design,
  method = c("within_block", "circular_shift", "global"),
  seed = NULL
)

# S3 method for class 'pair_rsa_design'
permute_labels(
  design,
  method = c("within_block", "circular_shift", "global"),
  seed = NULL
)
```

## Arguments

- design:

  An `mvpa_design`, `rsa_design`, or `pair_rsa_design` object.

- method:

  Character. One of `"within_block"` (default), `"circular_shift"`, or
  `"global"`.

- seed:

  Optional integer seed; RNG state is restored on exit.

## Value

A modified copy of `design` with the same class.

## Details

The generic dispatches on the design class:

- `mvpa_design`:

  Shuffles trial labels: `y_train`, `cv_labels`, and `targets` are
  replaced by permuted versions. `block_var` is never permuted.

- `rsa_design`:

  Permutes *item* labels, not RDM entries. The entries of an RDM are not
  independent observations: every item takes part in `n - 1` pairs, so
  shuffling the vectorised RDM would destroy that dependence and give an
  anti-conservative null. Relabelling items keeps it intact. The
  permutation is stored on the design as `item_perm`;
  `train_model.rsa_model` applies it to the rows of the neural pattern
  matrix before the brain RDM is computed. That is equivalent to
  permuting the rows and columns of every model RDM by the inverse
  permutation, so the model matrix and any cached kernel stay untouched.

- `pair_rsa_design`:

  As `rsa_design` in within-domain mode. For `pairs = "between"` the two
  item sets are permuted independently, each within its own blocks
  (`block_var_a`, `block_var_b`), and the second permutation is stored
  as `item_perm_b`.

When an RSA design excludes within-block pairs (the default whenever a
`block_var` is supplied and some pairs fall within a block), items are
exchangeable only within a block, and `method = "global"` is refused:
the null would compare neural patterns from the same run while the
observed statistic never does. If every block holds a single item, no
within-block shuffle can move anything and the method errors; when no
pairs are excluded in that situation, `method = "global"` is the valid
choice.
