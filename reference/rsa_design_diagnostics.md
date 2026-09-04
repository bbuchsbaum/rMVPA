# Diagnose the effective support of an RSA design

Summarises how much information an RSA design can bring to bear on each
predictor RDM. RDM entries are not independent observations: every item
enters `n_items - 1` pairs, so the `n_pairs` entries of a vectorised RDM
carry on the order of `n_items` independent pieces of information rather
than `n_pairs`. A regression coefficient's variance is inflated further
by the predictor's collinearity with the others, its variance inflation
factor (VIF). The ratio `n_items / VIF` is reported per predictor as the
effective number of items supporting that predictor's unique
contribution. Values below 10 mark a predictor whose coefficient will
vary across ROIs largely through noise. This is a heuristic screen, not
a test; it flags designs that cannot separate their predictors, it does
not certify those that can.

## Usage

``` r
rsa_design_diagnostics(design)
```

## Arguments

- design:

  An `rsa_design` or `pair_rsa_design`.

## Value

An object of class `rsa_design_diagnostics`: a list with `n_items`,
`n_pairs` (after any within-block exclusion), `n_predictors`,
`predictors`, `predictor_roles`, `max_abs_cor` and `max_abs_cor_pair`
(the most correlated pair of predictors), `vif` (named, `Inf` when the
design matrix is singular, `NA` for a constant predictor),
`items_per_predictor` (`n_items / vif`), and `threshold` (10).

## Details

For
[`pair_rsa_design`](https://bbuchsbaum.github.io/rMVPA/reference/pair_rsa_design.md)
objects in between-domain mode the item count is `n_a + n_b`, since each
item of either set enters every pair with the other set.

## See also

[`rsa_model`](https://bbuchsbaum.github.io/rMVPA/reference/rsa_model.md),
which computes this at construction and warns when a model predictor
falls below the threshold, and
[`run_permutation_searchlight`](https://bbuchsbaum.github.io/rMVPA/reference/run_permutation_searchlight.md)
for inference that respects the item-level dependence.

## Examples

``` r
set.seed(1)
items <- 12
shared <- matrix(rnorm(items * 4), items, 4)
a <- dist(shared)
b <- dist(shared + matrix(rnorm(items * 4, sd = 0.3), items, 4))
des <- rsa_design(~ a + b, list(a = a, b = b))
rsa_design_diagnostics(des)
#> RSA design diagnostics
#>   items: 12   pairs: 66   predictors: 2
#>   max |r| between predictors: 0.880 (a, b)
#>   effective items per predictor (items / VIF):
#>     a                    VIF    4.44   ~   2.7 items  <- below threshold
#>     b                    VIF    4.44   ~   2.7 items  <- below threshold
```
