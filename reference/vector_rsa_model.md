# Create a vectorized RSA model

This function integrates a vector_rsa_design and an mvpa_dataset to
create a vectorized RSA model.

## Usage

``` r
vector_rsa_model(
  dataset,
  design,
  distfun = cordist(),
  rsa_simfun = c("pearson", "spearman"),
  nperm = 0,
  save_distributions = FALSE,
  return_predictions = FALSE
)
```

## Arguments

- dataset:

  An `mvpa_dataset` object.

- design:

  A `vector_rsa_design` object.

- distfun:

  A `distfun` (distance function) for computing pairwise dissimilarities
  among image rows.

- rsa_simfun:

  A character string specifying the similarity function to use for RSA,
  one of `"pearson"` or `"spearman"`.

- nperm:

  Integer, number of permutations for statistical testing (default: 0).

- save_distributions:

  Logical, whether to save full permutation distributions (default:
  FALSE).

- return_predictions:

  Logical, whether to return per-observation similarity scores (default:
  FALSE).

## Value

A `vector_rsa_model` object (S3 class) containing references to the
dataset, design, and function parameters.

## Details

If \`return_predictions\` is TRUE, the output of \`run_regional\` or
\`run_searchlight\` will include a \`prediction_table\` tibble
containing the observation-level RSA scores.

## Examples

``` r
# \donttest{
  D <- as.matrix(dist(matrix(rnorm(5*3), 5, 3)))
  rownames(D) <- colnames(D) <- letters[1:5]
  labels <- factor(rep(letters[1:5], 2))
  block_var <- rep(1:2, each=5)
  des <- vector_rsa_design(D, labels, block_var)
  # mdl <- vector_rsa_model(dataset, des)
# }
```
