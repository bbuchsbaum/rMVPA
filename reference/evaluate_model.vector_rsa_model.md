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
if (FALSE) { # \dontrun{
  # Internal S3 method called during processing
  # perf <- evaluate_model(vector_rsa_model, observed, roi_data)
} # }
```
