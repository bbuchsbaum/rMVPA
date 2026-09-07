# Control parameters for the pattern model estimator

Control parameters for the pattern model estimator

## Usage

``` r
pattern_control(
  max_rank = 8L,
  x_scale = c("none", "sd"),
  y_scale = c("none", "sd"),
  noise = list(type = "diag_lowrank", rank = "auto", max_rank = 10L, shrink = 0.1),
  lambda_2 = 0,
  max_outer = 50L,
  tol = 1e-08,
  refine_path = FALSE
)
```

## Arguments

- max_rank:

  Maximum rank considered (capped at the eligible rank: number of
  classes minus one for categorical targets, the effective number of
  target dimensions otherwise, and never above the number of features or
  observations).

- x_scale:

  Feature scaling on training rows: `"none"` (centre only) or `"sd"`.

- y_scale:

  Continuous-target scaling before whitening: `"none"` or `"sd"`.

- noise:

  Residual covariance specification: a list with `type`
  (`"diag_lowrank"`, `"diag"`, or `"identity"`), `rank` (`"auto"`
  selects components above the Marchenko-Pastur edge, or an integer),
  `max_rank`, and `shrink` (shrinkage of the diagonal toward its
  median).

- lambda_2:

  Reserved for the penalized solver; must be 0. A ridge on the patterns
  turns the A-step into a Sylvester equation under a general residual
  covariance, so it is rejected rather than solved approximately.

- max_outer:

  Maximum number of alternating (C-step, A-step) updates.

- tol:

  Relative objective change that declares convergence.

- refine_path:

  Logical; when fitting a rank path, run the alternating refinement for
  every rank (default `FALSE`: path solutions are the exact reduced-rank
  optima, which coincide with the refined solution for the unpenalized
  objective).

## Value

A list of class `pattern_control`.

## Examples

``` r
pattern_control(max_rank = 3)
#> $max_rank
#> [1] 3
#> 
#> $x_scale
#> [1] "none"
#> 
#> $y_scale
#> [1] "none"
#> 
#> $noise
#> $noise$type
#> [1] "diag_lowrank"
#> 
#> $noise$rank
#> [1] "auto"
#> 
#> $noise$max_rank
#> [1] 10
#> 
#> $noise$shrink
#> [1] 0.1
#> 
#> 
#> $lambda_2
#> [1] 0
#> 
#> $max_outer
#> [1] 50
#> 
#> $tol
#> [1] 1e-08
#> 
#> $refine_path
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "pattern_control" "list"           
```
