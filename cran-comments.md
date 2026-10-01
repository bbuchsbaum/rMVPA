## Submission

This is a new submission.

rMVPA depends on neuroim2 (>= 0.20.0) and suggests fmridesign, fmrilss and
neurosurf. These packages, by the same maintainer, are submitted to CRAN
first; rMVPA is submitted once they are available, so no `Remotes` field
remains in the submitted DESCRIPTION.

## R CMD check results

0 errors | 0 warnings | n notes

* New submission.

* Imports includes 23 non-default packages. Each is used by the core
  analysis pipeline (dataset/design construction, cross-validation,
  searchlight iteration, result assembly); model families with optional
  back ends (e.g. glmnet, sda, pls, randomForest, e1071, xgboost) are in
  Suggests and are used conditionally via `requireNamespace()`.

* Installed package size (if reported): the bulk is R code (a large
  analysis API) and small example data in `inst/extdata` (a reduced Haxby
  et al. 2001 subject and the Kriegeskorte 92-image RDMs) used by the
  examples and vignettes.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.5.1
