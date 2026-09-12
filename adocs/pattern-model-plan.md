# pattern_model: active system and implementation plan

Status: existing estimator plus proposed agent-state/interoperability work.  
Date: 2026-09-12. Inspected base: `7c93f8d3afc9cafd4fcefaab22489de783601990` (`master`).

The full earlier assessment, implementation phases, deviations and dated receipts
are preserved byte-for-byte in
[the historical plan](archive/pattern-model-plan-2026-09-09.md). That document
records history, not a fresh instruction to rebuild completed phases. Old branch,
PR, timing and test-status statements must be interpreted at their recorded date.
This revision changes documentation, not numerical code or validation status.

## 1. Settled decisions and source of authority

The public name is **pattern_model**. Keep its pure numerical core inside rMVPA.
Retain the existing analysis-family integration through `create_model_spec`,
`fit_roi`, `roi_result` and scalar `output_schema`, with the thin global method.
Folds come from `model$crossval`, not `context$cv_spec`. The target accessor
completes the earlier cv_labels/targets architecture; do not change `y_train` to
return matrix targets or add another design hierarchy.

The joint control and evidence design is
[crossform's scientific-state-v1 contract](https://github.com/bbuchsbaum/crossform/blob/docs/agent-coherent-pattern-system/design/agent-system-contract.md).
Cross-package work is owned by
[the interoperability plan](https://github.com/bbuchsbaum/crossform/blob/docs/agent-coherent-pattern-system/design/pattern-model-interop-plan.md).
Its schema/examples are proposed interfaces, not callable R APIs. Existing code,
exported help, numerical contracts and source-bound tests govern behavior today.
Conflicts are reported and resolved through a dated decision, never by silently
choosing the newest prose. `pattern-agent-manifest.json` indexes this revision;
it does not duplicate a runtime capability registry.

The previously used `feature/sparse-smooth-pattern-model` ref did not resolve
through the connected repository during this revision. Changes are based on the
inspected master commit on a separate documentation branch, not a fabricated
working branch and not a direct modification of master.

## 2. Current state, not the old wishlist

| Area | Inspected/recorded state | Important boundary |
|---|---|---|
| Foundations | model_targets, matrix-free Haufe and spatial_graph are recorded as built | CV labels and model targets remain distinct; no dense p-by-p fallback |
| Core and spatial estimation | pattern_model/core/noise/prediction/resampling code is present | Whitening, dimensionless penalty definitions and monotone solver diagnostics follow the as-built record |
| Interpretation and locality | maps, immutable orthogonal views, regional predictions and stability are recorded as built | Local restriction is conditional on the shared task subspace, not an optimal independent ROI model |
| Confirmation | pattern_confirm, confirmation_plan, pattern_basis and component procedures are recorded as built | Error model and multiplicity scope are explicit; CR1 is approximate |
| Group loadings | R/pattern_group.R is present | Equal target subspaces and explicit one-to-one spatial identities; no arbitrary interpolation of uncertainty |
| Observation weights | Recorded as implemented on 2026-09-09 | Training/tuning weighted; reported prediction metrics unweighted; confirmation rejects nonuniform weights |
| Agent state and crossform export | Proposed by this revision | No new runtime API or capability is advertised as available |
| Supported rank, support envelope, local precision, joint hierarchy, response-set weighting | Deferred | No promotion from CV rank, favorable simulations or parameter placeholders |

Historical receipts are useful evidence pointers, but were not rerun here. The
active plan deliberately does not repeat old total test counts as a claim about
this new commit. Source presence, executed tests, hosted CI, merge and release
are separate states.

## 3. One central fitted relation

Keep the learned relation X = Y C A' + E as the numerical center. Retain its raw
and working coordinate transforms, noise model, rank-sized prediction core,
feature identities and training provenance. Prediction, model-implied geometry,
confirmed loadings and group evidence are different readouts, not interchangeable
objects with a generic importance field.

`remap_rrr_model` and `repmap_model` provide reduced-rank/cross-domain precedents;
the banded-ridge family and `feature_rsa_model` provide encoding/feature-analysis
and matrix-target precedents; `spacenet_tvl1` provides spatial decoding penalties.
Reuse their appropriate plumbing and benchmark them at matched estimands. Do not
route this model through scalar-target result wrapping or the conventional
classifier's duplicate-column screening. Do not claim ordinary PLS is a tested
thresholded-PLS comparator without a specified thresholded implementation.

For a reference task matrix T and neural measurement K, the induced geometry is
T(A'KA)T'. Export it as an external frozen prediction through a public crossform
contract, not by constructing `effect_geometry_fit` private fields. For Gaussian
prediction the same small core uses K=Psi^-1; downstream losses and inference
remain different.

## 4. Numerical and semantic invariants to preserve

Targets are centered/whitened on eligible training rows and mapped back through
stored transforms. The C-step's Procrustes argument depends on whitening; do not
reintroduce the earlier unwhitened derivation. Residual covariance is fixed from
a declared training pilot for an optimization stage. No silent per-iteration
covariance updates under an unchanged objective label.

The spatial A-step is a convex subproblem; the penalized factorization as a whole
is not thereby globally convex. Preserve convergence and restart diagnostics.
Use current dimensionless group-sparsity and graph-Laplacian signed-smoothing
scales. `support_smooth` is a different, deferred role, not a synonym for signed
quadratic smoothing. Warm starts across ranks were removed after convergence
problems; do not add them back as an untested caching optimization.

Vertex j is matrix column j, including basis-channel identity. Retain training-
only filtering and distinguish screened, outside-domain and estimated-zero
features. Identical values at two anatomical locations do not justify deleting
one column. Maps retain original units and transformations.

Keep fold-resolved and pooled ledgers. The latter gives one sorted observation
record; repeated assessment does not create independent observations. Report
predictive R2 against fold-training means, not squared correlation. Record the
actual target/observation/scoring weights separately.

Forward pattern A, calibrated-score Haufe A_H, decoding weights, signal SD and
conditional information have explicit names. The exact Haufe equality assumes
uncorrelated task and residual noise and the specified covariance. Singular score
spaces require the identifiable-subspace comparison. The Gaussian conditional-
information map is not empirical information for categorical labels or causality.

Local-restricted prediction uses (Psi_RR)^-1 and region-only test measurements.
Local-adapted covariance and independent ROI refits are different estimands.
An all-screened ROI yields the baseline; lack of retained fits yields a refusal.
Observed finite-sample accuracy is not guaranteed to increase with region size.

Confirmed loadings are unpenalized independent estimates in a frozen task basis,
not the sparse training map with added p-values. Current group loading inference
requires exact compatible target subspaces and available covariance. No arbitrary
Procrustes matching or many-to-one spatial transport of SEs.

## 5. Agent-first increments, not another seven-phase rebuild

### P-A0: accurate orientation and metadata-only inspection

Own the state adapters around existing validators/print methods, predominantly
`R/pattern_result.R`, `R/pattern_model.R` and `R/validate_analysis.R`. Derive human
and machine summaries from the same manifest. The brief names quantity, units,
source/fit/plan identity, usable readouts, missing prerequisites and bounded next
actions. Mark incomplete legacy provenance as unknown, never independent.

Exit: no neural data reads/fits during inspection; current call signatures are
used; unknown and proposed capabilities cannot become runnable; an agent with no
chat history finds the correct owning source and deferred feature status.

### P-A1: provenance and scientific plan differences

Extend, do not replace, existing observation IDs, basis IDs, fold hashes and
confirmation checks. Track transitive preprocessing/training ancestry and outcome
exposure. Distinguish a new scientific question from execution-only changes.
Rows renamed train/test do not establish independence. A fitted preprocessing
operator carries the rows that trained it.

Exit: source-alias, shared-normalization and inspected-test-set fixtures refuse
independent confirmation. Quantities with different weighting, conditioning,
metrics, baselines or correction families cannot be pooled silently. Do not claim
these checks establish physical independence.

### P-A2: bounded planning and retained execution

Expose fit counts, selected configuration, dominant copies, retained artifacts,
available local queries and a resource estimate with uncertainty. Inspect existing
fit paths before promising speedups. Add cancellation/atomic receipt and scoped
resume only for a demonstrated repeated workflow; preserve RNG and outcome state.
A changed plan forks rather than mutates a checkpoint. No global fit cache or new
scheduler is a prerequisite for using the model.

Exit: a repeated compatible geometry query does not refit; incompatible cache
keys refuse; partial work remains labeled partial; no dense graph/covariance
allocation or hidden whole-domain read in a regional query.

### P-A3: optional crossform interoperability

Implement the I0/I1/I2 sequence in the owning crossform interop plan: state views,
public frozen-prediction protocol, then retained partition-level confirmation
export where error capabilities match. Keep the adapter optional. A pooled
confirmation cannot invent independent partition estimates or residual sources.
Admit identity/bounded fixed metrics first; structured operator support is a
separate tested extension. Preserve current confirmation and group APIs.

### P-A4: cumulative knowledge and agent evaluation

Use the benchmark decision map below and existing history as initial experience
records. Bind every recommendation to a regime, receipt and invalidation trigger.
Run the crossform agent task battery with rMVPA entry points. Promote useful
failures into regression tests and reviewed recipes, not unconditional defaults.

Each increment needs one reviewer-readable diff, owning tests, independent
oracles where relevant, negative/refusal fixtures and a source-pinned receipt.
No calendar-duration promise substitutes for exit criteria.

## 6. Implementation boundaries and verification

For numerical changes, use an isolated R library and inspect exact available
filenames before running commands. Typical local gates are:

```sh
git status --short
git rev-parse HEAD
Rscript -e 'testthat::test_local(filter="pattern|global_analysis|fit_roi|plugin_extension_api|output_schema|collate")'
R CMD build .
# Then check the exact generated archive in an isolated library:
# R CMD check --no-manual <exact-archive-name>
```

Report exit status, warnings, skips and unrun gates. Do not overwrite the user's
library, regenerate all documentation for a prose-only change or call a historical
receipt current. Numerical/runtime tests are explicitly unrun in this revision.
Root Markdown and `adocs/` are already excluded from package builds.

For this revision, validate document links, JSON examples/schema, immutable archive
blob identity and changed-file scope. Link to the historical plan and benchmark
record rather than deleting them. Existing statistical test thresholds and code
are not changed.

## 7. Deferred decisions with evidence triggers

Local precision: revisit when a prespecified real-data comparison shows the
restricted head losing to a same-subspace local-covariance alternative. Separate
that from a missed-signal-subspace failure. TV support envelopes: revisit when a
specified localization objective justifies the complexity and signflip tests
remain acceptable. Supported rank: requires a validated higher-rank null, not a
new label for `rank_mean`. Joint hierarchy and target-set weighting require a
concrete use case and their own error/metric contracts.

The name, package placement and initial scope are settled. The next design work
is making existing capabilities legible and composable, not asking for those same
three decisions again.
