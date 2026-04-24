# Decorrelation of correlated model RDMs

Estimates layer- or model-specific RDM score vectors by removing shared
representational structure across a set of model RDMs. This is intended
as a preprocessing step before RSA when several model RDMs are strongly
correlated and the goal is to evaluate their more distinct
representational components.

## Usage

``` r
rdm_decorrelate(
  rdms,
  method = c("ordered_innovation", "zca"),
  similarity = c("spearman", "pearson"),
  objective = c("constraint", "penalty", "full_residual"),
  epsilon = 0.05,
  gamma_grid = seq(0, 1, by = 0.05),
  mu = 1,
  target = c("mean_abs", "mean_squared"),
  return = c("dist", "matrix", "vector"),
  project_correlation = FALSE,
  zca_ridge = 1e-08,
  max_combinations = 2e+05,
  tol = 1e-10
)
```

## Arguments

- rdms:

  Named list of square symmetric matrices or `dist` objects. All RDMs
  must have the same size and contain finite lower-triangular values.

- method:

  Character string, either `"ordered_innovation"` or `"zca"`.

- similarity:

  Character string. `"spearman"` rank-standardizes each RDM vector
  before decorrelation. `"pearson"` uses centered/scaled numeric RDM
  values directly.

- objective:

  Character string. `"constraint"` maximizes preservation subject to the
  requested decorrelation tolerance when feasible. `"penalty"` minimizes
  preservation loss plus `mu` times the decorrelation penalty.
  `"full_residual"` uses the fully residualized ordered innovations and
  ignores `epsilon`.

- epsilon:

  Non-negative tolerance for the mean pairwise adjusted-RDM association
  under `objective = "constraint"`.

- gamma_grid:

  Numeric vector in `[0, 1]` used for ordered innovation shrinkage. The
  first RDM is fixed at `gamma = 1`.

- mu:

  Non-negative penalty multiplier used when `objective = "penalty"`.

- target:

  Character string. `"mean_abs"` penalizes mean absolute pairwise
  correlation; `"mean_squared"` penalizes mean squared pairwise
  correlation.

- return:

  Character string. `"dist"` returns adjusted RDMs as `dist` objects,
  `"matrix"` returns symmetric matrices, and `"vector"` returns
  lower-triangular vectors.

- project_correlation:

  Logical. If `TRUE`, each adjusted RDM-score matrix is interpreted as
  `1 - correlation`, projected to the nearest positive semidefinite
  correlation matrix with unit diagonal using
  [`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html), and
  converted back to an RDM. This is off by default because projection
  can reintroduce RDM correlations.

- zca_ridge:

  Numeric ridge added only when the RDM-vector covariance is rank
  deficient under `method = "zca"`.

- max_combinations:

  Maximum number of exhaustive `gamma_grid` combinations for ordered
  innovation. Larger problems use coordinate search.

- tol:

  Numerical tolerance for degeneracy and rank checks.

## Value

An object of class `rdm_decorrelation_result` with adjusted RDMs,
lower-triangular vectors, selected shrinkage parameters, and
diagnostics.

## Details

This function works in lower-triangular RDM-vector space. With
`method = "ordered_innovation"`, RDMs are processed in the order
supplied: each RDM vector is decomposed into the component predicted by
earlier innovation vectors and a residual innovation. A shrinkage
parameter `gamma` controls how much of the predicted/shared component is
retained.

The default objective is constrained preservation: preserve as much of
each original standardized RDM vector as possible while forcing the mean
pairwise association among adjusted RDM vectors below `epsilon` when
feasible.

With `method = "zca"`, the standardized RDM vectors are symmetrically
whitened. This gives exact joint decorrelation in the Pearson sense when
the RDM-vector covariance matrix is full rank, but each output RDM is a
linear mixture of the original RDMs.

Ordinary identity shrinkage of a correlation matrix does not change
rank-based RSA correlations among off-diagonal RDM vectors, except at
the degenerate endpoint. This function instead shrinks shared RDM-vector
components.

## Examples

``` r
set.seed(1)
X <- matrix(rnorm(12 * 4), 12, 4)
r1 <- dist(X)
r2 <- dist(X + matrix(rnorm(12 * 4, sd = 0.2), 12, 4))
dec <- rdm_decorrelate(list(layer1 = r1, layer2 = r2), epsilon = 0.05)
dec$mean_abs_cor_after
#> [1] 8.121685e-16
```
