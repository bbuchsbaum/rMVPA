# Feature Selection Methods

Feature Selection Methods

## Methods

Two feature selection methods are available:

- FTest:

  One-way ANOVA F-test for each feature

- catscore:

  Correlation-adjusted t-scores using sda.ranking

## Cutoff Types

Two types of cutoffs are supported:

- top_k/topk:

  Select top k features

- top_p/topp:

  Select top p percent of features (0 \< p \<= 1)
