# Evaluate model performance for vector RSA

Computes the mean second-order similarity score and handles permutation
testing.

## Usage

``` r
evaluate_model.vector_rsa_model(
  object,
  predicted,
  observed,
  roi_data_for_perm = NULL,
  nperm = 0,
  save_distributions = FALSE,
  ...
)
```

## Arguments

- object:

  The vector RSA model specification.

- predicted:

  Ignored (vector RSA doesn't predict in the typical sense).

- observed:

  The computed second-order similarity scores (vector from train_model).

- roi_data_for_perm:

  New parameter

- nperm:

  Number of permutations from the model spec.

- save_distributions:

  Logical, whether to save full permutation distributions.

- ...:

  Additional arguments.

## Value

A list containing the mean RSA score (\`rsa_score\`), raw scores, and
optional permutation results (\`p_values\`, \`z_scores\`,
\`permutation_distributions\`).

## Examples

``` r
# Normally called internally during processing; shown here on toy scores.
ds <- gen_sample_dataset(c(4, 4, 4), nobs = 10, blocks = 2)
D <- as.matrix(dist(matrix(rnorm(5 * 3), 5, 3)))
rownames(D) <- colnames(D) <- letters[1:5]
des <- vector_rsa_design(D, factor(rep(letters[1:5], 2)), rep(1:2, each = 5))
mdl <- vector_rsa_model(ds$dataset, des)
trial_scores <- rnorm(10)
evaluate_model.vector_rsa_model(mdl, predicted = NULL, observed = trial_scores)
#> $rsa_score
#> [1] -0.1655465
#> 
#> $permutation_results
#> NULL
#> 
```
