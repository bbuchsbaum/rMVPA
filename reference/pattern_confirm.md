# Confirm a frozen pattern basis on independent observations

Confirm a frozen pattern basis on independent observations

## Usage

``` r
pattern_confirm(
  fit,
  dataset,
  design,
  block_var = NULL,
  inference = confirmation_plan(),
  observation_ids,
  discovery_ids,
  nuisance = NULL,
  feature_ids,
  preprocessing_id,
  subject_id = NULL
)

# S3 method for class 'pattern_confirmation'
print(x, ...)
```

## Arguments

- fit:

  A pattern fit, view, or global result with a refit.

- dataset:

  Confirmation observations by all input features, or an MVPA dataset
  whose training partition contains only confirmation observations.

- design:

  Confirmation targets (factor or numeric matrix), or an MVPA design
  with those targets in its training partition.

- block_var:

  Row-aligned independent run/block labels. Required for block-based
  error models; these labels do not add nuisance intercepts.

- inference:

  A
  [`confirmation_plan()`](https://bbuchsbaum.github.io/rMVPA/reference/confirmation_plan.md).

- observation_ids:

  Globally meaningful, unique confirmation row IDs.

- discovery_ids:

  Complete discovery row IDs in the same namespace, including all rows
  used for selection, tuning, or preprocessing. Required: the legacy
  fit's positional train/test IDs are not adequate provenance.

- nuisance:

  Numeric nuisance columns, without an intercept. Categorical nuisance
  variables must be coded by the caller using model.matrix.

- feature_ids:

  Unique IDs for all input columns, in their input order.

- preprocessing_id:

  A stable identifier for measurement units and the preprocessing
  recipe; matching strings are a caller assertion, not proof.

- subject_id:

  Optional unique participant ID, required for group analysis.

- x:

  A confirmation result to print.

- ...:

  Reserved.

## Value

A `pattern_confirmation`: original-feature-unit unpenalized estimates
and SEs, t and omnibus F tests, Holm-adjusted p values, complete
within-feature coefficient covariance (factorized for independent
errors), frozen basis and provenance. Sign-flip plans also return
bootstrap p and max-statistic adjusted p values. Nonestimable sampling
distributions are NA.

## Details

Regresses every input feature on frozen target scores and nuisance
columns, including an intercept. Discovery feature screening is not
repeated or used to select the confirmation hypothesis family. All
target scores must be estimable after nuisance adjustment. Independent
Gaussian homoskedastic errors give exact t/F tests; block methods are
approximate. Holm correction covers all feature-component tests, and
separately all omnibus tests. Component columns depend on the frozen
basis; omnibus tests are invariant to nonsingular changes of score
coordinates. Rank selection is not a rank hypothesis test and this API
does not report a supported population rank. IDs detect overlap but the
caller must establish genuine independence, including no shared
preprocessing fit or correlated repeated measurements.

## See also

`pattern_component_tests`, `pattern_basis`
