# Combine independent subject confirmations in shared target coordinates

Combine independent subject confirmations in shared target coordinates

## Usage

``` r
pattern_group(
  subject_confirmations,
  reference_basis,
  spatial_mapping = NULL,
  effects = c("random", "fixed")
)

# S3 method for class 'pattern_group_result'
print(x, ...)
```

## Arguments

- subject_confirmations:

  List of `pattern_confirmation` objects, each with a distinct
  subject_id and independent confirmation observations.

- reference_basis:

  A `pattern_basis` or discovery fit defining the shared raw target
  coordinates. Must be fixed independently of confirmation outcomes.
  Every subject basis must span exactly this target subspace.

- spatial_mapping:

  Optional list (one element per subject) of named character vectors:
  names are shared feature IDs, values are that subject's source feature
  IDs. Mappings must be one-to-one and cover the same shared features.
  With NULL, all subjects must have the same feature-ID set. This is
  correspondence, not spatial interpolation or parcel averaging.

- effects:

  `"random"` estimates the equally weighted population mean across
  subjects; `"fixed"` estimates a common effect using full
  inverse-covariance weighting.

- x:

  A group result to print.

- ...:

  Reserved.

## Value

A `pattern_group_result` with mean estimates, full mean covariance,
component tests, invariant omnibus tests and effect norms, moment
estimates of between-subject covariance, aligned subject estimates,
leave-one-subject-out descriptive expression, prediction summaries, and
provenance. Save the complete object with saveRDS.

## Details

Subject coefficients and their full covariance are transformed to the
reference basis before pooling. Unequal target subspaces are rejected;
Procrustes approximation would change the estimand. Original feature
units, target names/units, nuisance meaning, and preprocessing recipes
must agree. The preprocessing identifier and nuisance column names are
checked; their scientific equivalence remains the caller's
responsibility. Target sampling and omitted target effects must also
permit a common conditional estimand.

Random effects use the arithmetic subject mean and sample coefficient
covariance divided by number of subjects. Between-subject covariance is
the positive-semidefinite part of sample covariance minus average
sampling covariance (a moment estimate, not REML). Sampling error
remains in the mean covariance; it is not added a second time. Component
t tests have s-1 df; the omnibus uses Hotelling's T-squared transformed
to F(r, s-r). These are exact for iid Gaussian subject estimates with a
common total covariance, and approximate with heterogeneous
within-subject precision or nonnormal effects. Requires s \> r+1.
Subject count, not observation count, sets df.

Fixed effects use multivariate GLS and asymptotic normal/chi-squared
tests treating estimated within-subject covariances as known;
heterogeneity Q has r(s-1) df under the common-effect null. They do not
support population generalization. Holm correction is separate for all
component and all omnibus tests. Singular mean covariance yields NA
omnibus inference.

Norms, omnibus tests, heterogeneity trace, and expression are invariant
to orthogonal reference rotations. Component estimates and tests are
basis dependent. Arbitrary scaling of reference axes changes norms.
Expression projects a subject's coefficient vector onto the normalized
mean of the other subjects; it is descriptive, not an independent group
prediction. Interpolation would need cross-feature covariance, which
confirmations do not store, so many-to-one spatial mappings are
deliberately rejected.

## See also

`pattern_confirm`, `pattern_basis`
