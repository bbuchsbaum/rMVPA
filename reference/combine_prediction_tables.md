# Combine prediction tables

Combines multiple prediction tables (e.g., from different models or
regions) into a single table. Supports weighted combination and
collapsing regions.

## Usage

``` r
combine_prediction_tables(
  predtabs,
  wts = rep(1, length(predtabs)),
  collapse_regions = FALSE
)
```

## Arguments

- predtabs:

  A list of prediction tables (data frames) to be combined.

- wts:

  A vector of weights, with the same length as `predtabs`. Default is
  equal weights.

- collapse_regions:

  A logical value; if TRUE, regions are collapsed in the final
  prediction table.

## Value

A combined prediction table (data frame).

## Details

For classification, this function pools class-probability columns
(\`prob\_\*\`) and returns pooled predicted class and correctness. For
regression (numeric \`observed\`), it pools \`predicted\` values and
returns residual diagnostics.

## Examples

``` r
# Create example prediction tables
observed = factor(sample(letters[1:2], 10, replace = TRUE))
predtab1 <- data.frame(.rownum = 1:10,
                       roinum = rep(1, 10),
                       observed = observed,
                       prob_A = runif(10),
                       prob_B = runif(10))
predtab2 <- data.frame(.rownum = 1:10,
                       roinum = rep(2, 10),
                       observed = observed,
                       prob_A = runif(10),
                       prob_B = runif(10))

# Combine the tables
combined_table <- combine_prediction_tables(list(predtab1, predtab2))
```
