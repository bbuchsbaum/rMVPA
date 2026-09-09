# Test frozen component association and incremental decoding value

Test frozen component association and incremental decoding value

## Usage

``` r
pattern_component_tests(fit, confirmation, dataset, design, calibration = NULL)
```

## Arguments

- fit:

  The discovery fit used for confirmation.

- confirmation:

  A `pattern_confirmation` from that fit.

- dataset, design:

  The exact confirmation data and targets used by `pattern_confirm`.
  Content hashes are checked before computing scores.

- calibration:

  Optional list with `dataset`, `design`, and `observation_ids`. Its IDs
  must be included in the declared discovery set. Full and
  leave-one-component-out ordinary least-squares decoding heads are
  fitted on these rows alone. Without calibration only association is
  tested.

## Value

A list with an association table and, when calibration is supplied, an
incremental-loss table, fitted decoding heads, and row-level losses.

## Details

Association is the partial correlation of each frozen decoded score with
its corresponding frozen target score, adjusting for confirmation
nuisance columns (not for other components). Tests use the confirmation
error model and Holm correction over components.

Incremental gain is reduced-head minus full-head squared prediction
error, averaged over raw target columns, then within independent blocks,
then equally across blocks. Categorical targets use one-hot squared
loss; linear heads are not probability models. Positive gain favors
retaining the component. The paired block-mean t test is approximate;
sign-flip plans use a centered, studentized block wild bootstrap, also
approximate. These tests concern separately fitted decoding heads, not
the original Gaussian posterior decoder or a supported latent rank.
Component tests depend on the frozen basis. Score degeneracy is reported
for association and rejected when a decoding head cannot be estimated.
