# Issue 92: joint endpoint decomposition validation

Date: 2026-11-04. Issue: <https://github.com/bbuchsbaum/rMVPA/issues/92>.
Base HEAD: `3dd1ac725963609e58331fcd9d9d1497f0e912f8`, on
`fix/cran-blockers`. Validation was performed on the uncommitted implementation;
the subsequent commit and push are recorded in Git history.
The pre-existing `tests/testthat/testthat-problems.rds` was preserved.

## Delivered behavior

`era_endpoint_decomposition()` fits complete, diagonal-excluded joint endpoint
effects with explicit one-to-one keys and competitor blocks. Templates are rows,
retrieval patterns columns. The default pattern-sum path avoids a quadratic
similarity allocation; explicit matrices are optional. Small blocks, missing or
constant patterns, and unsupported weighted designs have explicit statuses or
refusals. No generalized missing-pair/weighted solver is claimed.

`summarize_era_endpoints()` retains the raw shared baseline, structural endpoint
zeros and descriptive RMS. Behavioral contributions use common complete cases,
one nuisance projection and one Pearson denominator. Separate Pearson or
Spearman correlations remain separately named and guard near-constant inputs.
Common-row counts, exclusion reasons and identity errors are retained.

`era_endpoint_spec()` opts ERA and feature RSA into scalar regional/searchlight
outputs. ERA additionally accepts `era_association_score = "joint_excess"`.
Feature scoring uses unchanged OOF predictions, records their fold identities,
and rejects repeated/partial assessments or competitor pools crossing folds.
Regional item retention is optional. Existing scores and tuning are unchanged.
The [methods article](../vignettes/articles/ERA_Endpoints.Rmd) describes the API,
geometry, interpretation and inference boundaries.

## Dependency correction

Commit `8d601883640fabebc9308250952d1cd167bc8e28` introduced the minimum
`neuroim2 >= 0.20.0` to avoid its masking `scale()` export. Current verification
found installed/CRAN 0.19.1, GitHub master 0.19.1.9000, and local
`~/code/neuroim2` 0.17.0 at `9aa1bcb81340a68ed7da7ef8856706a8fc359d7d`.
Installed 0.19.1 does not export `scale`; the local commit also explicitly removes
that generic. The dependency floor is corrected to tested 0.19.1. No scientific
tolerance or existing evaluation threshold was relaxed. The old CRAN-check
receipt is explicitly historical.

## Observed validation

R 4.3.3 on arm64 macOS, with an isolated package library at
`/private/tmp/rmvpa-issue92/lib`. Existing dependency libraries were read only.
The user's R library and neuroim2 checkout were not changed. Temporary R/Rscript
launchers ensure child workers use R 4.3 rather than the system's R 4.6; a
two-worker smoke check verified R version, library path and neuroim2 version.

- Scoped source tests: **122 tests, 1,400 passing expectations, zero failures,
  zero errors, zero skips, 22 warnings**. Warnings comprise 14 existing ERA
  metadata warnings, one existing constant-prediction warning, and seven Git
  metadata warnings from exporting in a non-Git temporary snapshot.
- Scope: `era_endpoint`, `era_association`, `era_rsa_model`, `era_rsa_design`,
  `era_rsa_searchlight_fast`, `feature_rsa_predictions`,
  `feature_rsa_fold_metrics`, `output_schema`, `collate`, `save_results`.
- The two new test files also passed against the **installed** package:
  **303 passing expectations**, zero failures, one temporary-snapshot Git warning.
- Independent off-diagonal `lm` oracle; residual row/column sums; additive-matrix,
  diagonal-injection and transpose controls; randomized matrix/sum equivalence;
  separate phase reorderings; unequal block sizes; constant retrieval/encoding;
  common-case nuisance projection; separate partial-rank correlation oracle;
  small/nonfinite/constant-input refusals; near-constant behavioral subsets.
- Identical voxel sets agree across regional, general-purpose searchlight and
  dedicated ERA engines. OOF feature predictions, observations, folds, tuning
  diagnostics and existing metrics are unchanged. NIfTI round trips preserve
  mean excess, structural NA correlations and local item counts.
- `R CMD build --no-build-vignettes`: exit 0.
- `R CMD INSTALL` into the isolated library: exit 0, including namespace and
  shared-library load checks.
- `R CMD check --no-manual --no-vignettes --no-tests` with
  `_R_CHECK_FORCE_SUGGESTS_=false`: exit 0; **zero errors, two warnings, two notes**.
  Both warnings concern the deliberately unbuilt seven existing package
  vignettes. Notes concern unavailable optional suggested packages and installed
  size (5.3 MB). Syntax, namespace, dependencies, code analysis, Rd, compiled
  code and examples passed.
- New article rendered successfully against the installed package. Pandoc
  emitted a deprecated `--highlight-style` option warning.
- `git diff --check`: passed.

Archive SHA-256:
`be1c647f743426033cc33add3ffd4becf459ce0c861e59c349d41d5697d8e332`.
[Source hashes](era-endpoints-issue92-source.json) identify the 19 changed/new
implementation, test, help and documentation files in the validated snapshot.
This receipt and the hash manifest are documentation added after validation.
Local logs, archive, rendered HTML and isolated library are under
`/private/tmp/rmvpa-issue92/`; these temporary artifacts may be removed later.

## Unrun and next action

The full package test suite, all existing vignette builds, manual, hosted CI,
real-data validation, group inference and endpoint permutations were not run.
No participant-derived data were used or published. Group inference and
permutation nulls are not implemented by this feature; existing feature
permutation outputs continue to cover only their original metrics.

Smallest remaining validation action: run the full release checks with the
intended optional dependencies before a CRAN submission. This receipt does not
claim hosted CI, PR review, issue closure or release qualification.
