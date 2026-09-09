# Evaluate regional access under a whole-brain pattern fit

Evaluate regional access under a whole-brain pattern fit

## Usage

``` r
local_performance(result, regions, independent_roi = NULL)
```

## Arguments

- result:

  A `pattern_global_result` retaining fold fits.

- regions:

  Named list of input-column positions or logical masks. These are
  dataset matrix columns, not voxel IDs. Regions may overlap.

- independent_roi:

  Optional named list of fold-resolved `pattern_ledger`s from
  independent ROI models on exactly the same rows and folds. Truth,
  target type, partition, and baseline must agree.

## Value

A table with region, metric, `whole_brain`, and `local_restricted`
columns (and `independent_roi` if supplied). Regional pooled and
fold-resolved ledgers are retained as attributes.

## Details

Each region uses the whole-brain task representation and the exact
marginal covariance `Psi_RR`; neither patterns nor noise are refitted.
All preprocessing is columnwise. Repeated assessments are averaged per
observation before scoring, just as for the whole-brain ledger.
R-squared uses fold training means. This measures locally accessible
information under the shared model, not optimal performance of a
separately fitted ROI. Finite-sample regional accuracy need not be below
whole-brain accuracy.
