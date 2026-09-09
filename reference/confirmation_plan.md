# Specify a confirmation error model

Specify a confirmation error model

## Usage

``` r
confirmation_plan(
  error = c("independent", "block_robust", "sign_flip"),
  n_resamples = 999L,
  seed = 1L
)
```

## Arguments

- error:

  Independent Gaussian rows, a CR1 sandwich over independent blocks, or
  a restricted-residual block wild bootstrap with Rademacher signs.

- n_resamples:

  Number of Monte Carlo sign draws (at least 99).

- seed:

  Nonnegative integer seed. Resampling restores the caller's RNG.

## Value

A `pattern_confirmation_plan`.

## Details

The independent model uses residual degrees of freedom. Block sandwich t
and Wald F tests use number of blocks minus one and are small-sample
approximations. The sign-flip option uses the same sandwich statistics
and refits residuals under each null. With estimated nuisance effects
this is an approximate wild bootstrap, not an exact permutation test.
Independent blocks and enough blocks for the tested dimension are
required. No option establishes independence from row labels alone.
