# pattern_model: evidence, resource decisions and agent benchmarks

Status: active decision map plus proposed new evaluations. Date: 2026-09-12.
The full numerical tables and dated receipts are preserved byte-for-byte in
[the historical benchmark record](archive/pattern-model-benchmarks-2026-09-09.md),
from rMVPA `7c93f8d`. Those experiments were **not rerun for this revision**.
The current implementation plan is [pattern-model-plan.md](pattern-model-plan.md).

## 1. Read an experiment as a bounded experience record

A reusable record names source/artifact revisions, generator/data, observation
and feature counts, target rank, blocks, signal organization, covariance regime,
comparison estimand, selection rules, hardware/BLAS, resource measurements,
uncertainty, conclusion, limitations and the condition that would invalidate its
recommendation. Preserve negative findings and competing explanations. One
simulation family cannot decide a universal default.

| Record | Recorded finding, not a fresh result | Action justified | Boundary / revisit trigger |
|---|---|---|---|
| PM-E01 | At p=100k, n=400, auto-sparsity workflow: 77.5 s and roughly 2.2 GB on the recorded single-core environment | Profile full-data copies/decompositions before optimizing a tiny prox kernel | Hardware, n/r, CV/path, penalty and retention changes require new measurements |
| PM-E02 | Penalized A-step about 3% of one measured fit | Do not add C merely because the model is whole-brain | New penalty or profile showing a materially different bottleneck |
| PM-E03 | Shared-pattern regional covariance heads close on the planted simulation | Keep local-precision extension conditional | Prespecified real-data/head-only comparison demonstrates a reproducible gap |
| PM-E04 | Dense weak signal: auto sparsity worse than unpenalized model (0.471 vs 0.738 in the recorded case) | Include an unpenalized comparator; treat sparsity as a scientific assumption | Additional scoped simulations/data may alter the recommendation, not erase this result |
| PM-E05 | Signed smoothing improved smooth-pattern recovery while worsening decoding in one scenario | Score localization and prediction separately | Different objective or signal organization |
| PM-E06 | CR1 nominal 5% omnibus rejection 8.4%; block wild bootstrap 4.1% in the recorded null experiments | Display CR1 as approximate and surface the tested limitation | Different block/sample/error regimes need calibration; bootstrap is not universally exact |
| PM-E07 | Group moment procedure near nominal in its independent Gaussian sufficient-statistic simulation | Use only the declared compatible-coordinate/group assumptions | Arbitrary heterogeneity, interpolation or missing covariance is outside that evidence |

PM-E IDs are local experience-record identifiers, not new crossform scientific
claim promotions. They link to the historical scripts and observations. Numerical
and inferential consequences must remain separated from agent-product performance.

## 2. Correct the comparison questions before running baselines

Whole-brain classifiers, independently fitted ROI classifiers, a selected single
searchlight and a map of searchlight scores are different outputs. Compare
prediction on identical assessment rows and frozen inner-selection rules; compare
localization on its own known-truth target. An oracle test-selected sphere is an
explicit retrospective diagnostic, never the headline fair competitor.

To diagnose local loss, hold A_R fixed while changing only the covariance head.
Then compare the best eligible head with an independently fitted local model.
Do not conclude covariance is harmless globally because a small simulation found
no difference. Plain PLS is not thresholded PLS. Record each baseline's actual
implementation, version, preprocessing, tuning space and failure handling.

## 3. Resource receipts for the entire workflow

Account separately for source reads, copies, target whitening, pilot covariance,
rank/penalty search, refit, global prediction, retained fold artifacts, regional
queries, bootstrap/group inference and saving. Report cold and compatible-reuse
runs, measured RSS versus buffer estimates, thread count and RNG. Include failure,
cancellation and resume cost. Partial results do not qualify as complete timings.

No cross-rank warm starts without new correctness evidence. No cache reuse across
training/label/transform changes without a proved equivalence. Profiling, not an
aspirational abstraction, decides whether persistent caching or compiled kernels
are worthwhile.

## 4. Agent-product acceptance

Use the shared
[agent evaluation plan](https://github.com/bbuchsbaum/crossform/blob/docs/agent-coherent-pattern-system/design/agent-evaluation-plan.md).
Run the same tasks with the old interface and the proposed metadata view. Record
correct completion, scientific mistakes, unnecessary fits/reads, tool calls,
context, wall time, peak memory and human corrections separately. Efficiency is
compared only among correct completions. No desired speedup is asserted in advance.

Mandatory traps include stale branches/receipts, matrix targets, weighted training
versus unweighted scoring, missing retained folds, cropped precision, rank-zero
regions, outcome-exposed confirmation, incompatible subject bases and forged
crossform fit classes. Each has an expected refusal or exact correct action.

## 5. Executable sources and promotion

Existing recorded numerical scripts are under `inst/benchmarks/pattern_model/`:
`bench_pattern_model.R`, `bench_vs_baselines.R`, `validate_confirmation.R`, and
`validate_group.R`. Verify their current contents and environment before rerunning.
Store new outputs outside tracked source unless intentionally publishing a reviewed
sanitized receipt. Do not install dependencies into the user's library by default.

A new run appends a scoped record and source fingerprint. It never overwrites the
old experiment or retroactively edits eligibility/thresholds to pass. An adaptive
exploration becomes a new prospective question only with eligible untouched data
and a genuinely prior protocol. Publishing schema-valid documentation does not
promote any numerical or scientific evidence status.
