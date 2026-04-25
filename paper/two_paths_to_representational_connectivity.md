---
title: "Two paths to representational connectivity in rMVPA"
author: "Bradley Buchsbaum"
date: "April 2026"
---

> **Draft methods note** accompanying the
> [`rMVPA`](https://github.com/bbuchsbaum/rMVPA) R package.
>
> The note formalises the pair-observation model-space RSA framework
> implemented in `rsa_model(..., return_fingerprint = TRUE)` and
> `model_space_connectivity()`, contrasts it with the feature-RSA
> connectivity path implemented in `feature_rsa_model()` and
> `feature_rsa_connectivity()`, and describes a memory-bounded extension
> for searchlight analyses based on k-means anchors.

## Abstract

Representational similarity analysis (RSA) asks how strongly a region's
neural pattern geometry tracks a model-derived dissimilarity structure.
Once one has fitted RSA in many regions, the natural follow-up is
*representational connectivity*: which regions share that geometry, and
how do you ask the question across two stimulus sets, or across whole-brain
searchlight, without materialising an `n_centers × n_centers` matrix?

We argue that within-unit RSA and across-unit representational
connectivity should be two views of the same fitted object. Concretely, a
single call to `rsa_model(..., return_fingerprint = TRUE)` produces, per
analysis unit (an ROI, a searchlight sphere, or a domain-specific bin), a
small whitened projection of the unit's neural pair-dissimilarity vector
onto an orthonormal basis of the model-RDM subspace — the unit's
*model-space fingerprint*. Stacked into a matrix `F`, the fingerprints
yield ROI-by-ROI connectivity as `F Fᵀ`, with no further fitting and at
memory cost O(K) per unit where K is the number of model RDMs.

We compare this *projection* path to the alternative *learned* path
implemented in `feature_rsa_connectivity()`, which fits a
cross-validated map between neural patterns and a feature space. The two
paths answer subtly different questions and have very different
computational footprints; choosing between them turns on whether you
have an explicit model RDM (projection) or a feature space whose
representation must be learned (learned). For whole-brain searchlight we
introduce a k-means anchor scheme that summarises representational
connectivity at memory O(n_centers × k), avoiding the prohibitive
`n_centers²` materialisation.

## 1. Background: RSA, briefly

Let `X ∈ ℝ^{T × V}` be a region's neural pattern matrix (T trials or
conditions × V voxels). Let `R = D(X) ∈ ℝ^{T × T}` be the
representational dissimilarity matrix induced by some distance metric
`D` (typically `1 − Pearson correlation` or correlation-distance).

Standard RSA compares the lower-triangular vector `r = vec_lt(R)` of
length `p = T(T−1)/2` with one or more model RDM vectors
`m₁, …, m_K`. Letting `M = [m₁ … m_K]`, RSA fits

$$ r ≈ M \beta + \varepsilon $$

and reports either each $\beta_k$, the partial correlation
$\text{cor}(r, m_k \mid m_{≠k})$, or just the marginal Spearman
correlation $\rho(r, m_k)$. Whatever the choice of estimator, RSA
returns a vector of *K scores per region* describing how much of that
region's geometry each model RDM explains.

When one has run RSA in `n` regions, the resulting `n × K` score matrix
already implicitly contains a representation of inter-region
relationships: regions that load on the same axes of the model space
should produce similar score vectors. *Representational connectivity*
makes that observation explicit and adds two refinements: (i) the score
axes should be orthogonalised so that comparing rows is well-defined,
and (ii) the comparison should be expressed at second order via an
inner product rather than via raw Pearson correlation across small
score vectors.

## 2. The pair-observation view

The shared abstraction is **pairwise model-space representation per
analysis unit**. Each "unit" `u` (an ROI, a searchlight, a parcel × subject
cell, …) supplies:

- a length-`p` neural pair vector `r_u = vec_lt(R_u)`, and
- the same model design matrix `M ∈ ℝ^{p × K}` shared across units.

For an ROI analysis the units are the ROIs; for a searchlight, the
units are searchlight centres; for a cross-subject design the units are
(subject × session) cells. The pair design `M` may itself be
non-trivially constructed: lower-triangular within-domain pairs are the
classical case, but rectangular `n_a × n_b` between-domain pair tables
(`pairs = "between"` in `pair_rsa_design()`) and function-valued model
entries (a function `f(item_i, item_j) ↦ ℝ` evaluated over the pair
table) are equally first-class. The remainder of this note is agnostic
to which pair design produced `M`.

## 3. Fingerprints: whitened projection onto the model-RDM basis

Let `M = QΣVᵀ` be the (thin) singular value decomposition of `M` over
the cells used in the analysis (after standardisation: column-centred and
unit-variance for Pearson, rank-then-standardised for Spearman). The
columns of `Q ∈ ℝ^{p × rank(M)}` are an orthonormal basis of the
column space of `M`.

The **fingerprint** of unit `u` is

$$ f_u = \frac{1}{\sqrt{p − 1}} \, Q^\top y_u, $$

where `y_u` is the standardised version of `r_u` under the same
transformation applied to `M`. Geometrically, `f_u` is the
coordinates of unit `u`'s neural pair vector in the orthonormal axes of
the model space. Each entry of `f_u` is the standardised RSA score on
one orthogonal axis — a familiar number, just expressed in a basis that
makes downstream algebra clean.

Three properties matter:

1. `‖f_u‖² = ‖P_M y_u‖² / (p − 1)`, the **model-explained strength**:
   how much of the unit's standardised pair vector lies in the model
   subspace.

2. The cosine of `f_u` against `f_v` is the *profile similarity*:
   how much the two units share the same orientation in model space,
   irrespective of strength.

3. `f_u` lives in `ℝ^K` (or lower-dimensional if `M` is rank-deficient).
   *K typically lies between 2 and 10*, even for whole-brain
   searchlight. This is the compression that makes downstream
   connectivity tractable.

In rMVPA, fingerprints are produced by `rsa_model(..., return_fingerprint
= TRUE)` and travel with the result of `run_regional()` /
`run_searchlight()` as `attr(result, "fingerprints")$scores`.

## 4. Representational connectivity = F Fᵀ

Stack the unit fingerprints into a matrix `F ∈ ℝ^{n × K}` (rows = units).
Then

$$ S = F F^\top \in ℝ^{n × n} $$

is the **representational connectivity matrix** through the model-RDM
subspace. The `(u, v)` entry is the inner product of the two units'
model-space coordinates: a strength-sensitive similarity that reduces
to the squared norm `‖f_u‖²` on the diagonal (model-explained strength)
and to a strength-weighted cosine off-diagonal.

Several decompositions of `S` carry different scientific meaning:

- **Profile similarity**: `Sᶜ = D F (D F)ᵀ` where `D = diag(‖f_u‖⁻¹)`.
  Cosine similarity of fingerprints; orientation only.
- **Component similarity**: `S = Σ_k F_k F_kᵀ`, where `F_k` is the
  k-th column of `F`. Each rank-1 matrix shows the connectivity
  contribution of one orthogonal model axis. Under `basis = "pca"` the
  first axis is typically the *common* axis (shared geometry across
  units); the rest are *differentiating* axes.
- **Raw vs residual similarity**: `Y Yᵀ / (p−1) = S + S_residual`, where
  the residual sits in the orthogonal complement of the model subspace.
  Useful for testing whether two units share *additional* geometry
  beyond your declared models.

`rdm_model_space_connectivity()` returns all of these as named matrices
on a single object, plus diagnostics (eigenvalues, axis-to-original-RDM
correlations).

## 5. Two paths, summarised

rMVPA exposes two complementary routes to representational connectivity.
The shared question is: *do regions A and B share representational
geometry?*

| | **Model-space connectivity** | **Feature-RSA connectivity** |
|:--|:--|:--|
| Implementation | `rsa_model(..., return_fingerprint = TRUE)` → `model_space_connectivity()` | `feature_rsa_model(..., return_rdm_vectors = TRUE)` → `feature_rsa_connectivity()` |
| What you supply | Explicit model RDMs | A feature matrix `F` (or similarity `S`) |
| What is fitted | Nothing; project onto a fixed basis | Cross-validated PLS/PCA/glmnet maps neural patterns to features |
| What is compared | Whitened projections in model space (length K) | Predicted trial-by-trial RDM vectors (length `p`) |
| Sensitivity | Restricted to the model subspace | Anything the feature space can capture |
| Memory per unit | O(K), typically 2–10 numbers | O(p), typically 10⁴–10⁵ numbers |
| Searchlight cost | Tractable via k-means anchors (§6) | Marginal; needs file-backed RDM storage |
| Best when | You have specific theoretical RDMs | You don't have an explicit RDM but do have features |
| Diagonal of `S` | Model-explained strength | Within-unit CV alignment |

A useful rule of thumb: **declared model → model-space connectivity;
learned model → feature-RSA connectivity.** The two are not redundant.
Where they disagree is informative: model-space connectivity inherits
the strengths and limitations of your declared models, while feature-RSA
connectivity inherits the inductive biases of its fitting procedure.

## 6. Searchlight: anchors via k-means

Whole-brain searchlight has on the order of `n ∼ 10⁵` centres. The
fingerprint matrix `F ∈ ℝ^{n × K}` itself is small (`10⁵ × 10` doubles ≈
8 MB), but `F Fᵀ ∈ ℝ^{n × n}` would be ~80 GB — not materialisable.

We replace `F Fᵀ` with a low-rank summary by selecting `k` *anchor*
searchlights and reporting an `n × k` similarity matrix between every
centre and each anchor. Letting `A ⊂ {1, …, n}` be the anchor index
set with `|A| = k`,

$$ \tilde S = F F_A^\top \in ℝ^{n × k}. $$

Each column of `\tilde S` becomes a brain map — voxel by voxel,
"how similar is this searchlight's representational geometry to the
geometry at anchor `j`?" — built directly via the existing
`build_output_map()` pipeline.

Anchor selection is the design choice. Three strategies, in order of
sophistication:

1. **k-means in fingerprint space** (default in
   `model_space_connectivity.searchlight_result`): cluster the rows of
   `F` into `k` clusters and pick the searchlight whose fingerprint is
   closest to each cluster centroid as the anchor. Each anchor map says
   "how similar is this voxel to one of `k` representative pattern
   types." Interpretable, parameter-light.

2. **Leverage / CUR sampling**: probability ∝ `‖f_u‖²` (or row leverage
   of `F`). Gives a Nyström approximation of `F Fᵀ` with provable
   spectral preservation. Theoretical optimum when one cares about
   reconstructing the full connectivity spectrum.

3. **Explicit seeds**: the user supplies anchor center IDs (e.g. peak
   coordinates from a prior analysis or specific parcels). Bypasses
   clustering entirely.

Memory cost is `O(n × k)` regardless of `K` or the number of trials.
For typical `k = 20` anchors and `n = 10⁵` centres this is 16 MB —
fits trivially. We default to k-means anchors with cosine
("`scale = "norm"`") similarity for interpretability, exposing CUR as
an option for users who want spectral guarantees, and explicit seeds
for users with prior knowledge.

The k-means path also has a graceful clamping: if fewer than `k`
distinct fingerprints are present, `k` is reduced to the number of
distinct rows; if `k` equals or exceeds `n`, the function skips
clustering and returns the full `n × n` similarity directly.

## 7. Storage of pair vectors at scale

A separate concern from the connectivity summary is *storing* the
neural pair vectors themselves. For 500-trial designs the lower
triangle has ≈ 125 000 cells. Across 10⁵ searchlights this is 100 GB —
again, not materialisable.

Three complementary strategies are available in rMVPA:

- **Fingerprints not RDMs.** The `return_fingerprint` path discussed
  above replaces each unit's 125 000-cell pair vector with a length-K
  fingerprint. Memory drops from 100 GB to ~8 MB. This is the
  recommended default.

- **Random projection (Johnson–Lindenstrauss) of pair vectors.** When
  RDM-shaped output is genuinely required (e.g. for cross-RSA, model
  selection downstream), each unit's pair vector can be sketched onto
  a fixed `m`-dimensional subspace via a sub-Gaussian sketch
  `S ∈ ℝ^{p × m}`. Pair-cosine and pair-correlation are preserved up
  to JL distortion at memory `O(n × m)`; for `m = 1000` and `p = 125000`
  this is a 125× compression.

- **File-backed batching.** For full-fidelity preservation (no
  compression), `feature_rsa_model(save_rdm_vectors_dir = …)` writes
  per-batch RDS files; `feature_rsa_connectivity()` consumes them
  block-wise. Trades RAM for I/O.

## 8. Implementation in rMVPA

The framework is implemented across a handful of functions:

| Function | Role |
|:--|:--|
| `pair_rsa_design()` | General pair-coordinate system: within (lower-tri), between (rectangular A×B), function-valued model entries, optional nuisance RDMs. Inherits from `rsa_design`. |
| `rsa_model(..., return_fingerprint = TRUE)` | Per-unit fit; on top of standard RSA scores, projects the neural pair vector onto an orthonormal basis of the model-RDM subspace and returns the standardised coordinates. |
| `train_model.rsa_model()` | Branches on `pair_kind` (within vs between) for neural-side pair extraction. |
| `model_space_connectivity()` | Unified entry point. S3 dispatch on `regional_mvpa_result` (fingerprint-driven), `searchlight_result` (k-means anchors), and a default that forwards to `rdm_model_space_connectivity()` for matrix inputs. |
| `rdm_model_space_connectivity()` | Matrix-level workhorse: standardise → orthonormal basis → project → connectivity. Accepts ROI RDM matrices, lists, or feature-RSA result tibbles. |
| `feature_rsa_connectivity()` | The learned-feature alternative. Operates on per-ROI predicted RDM vectors; reports cross-ROI Spearman/Pearson similarity and asymmetric source→target generalisation. |

Two real-data vignettes demonstrate the framework on canonical public
data: `vignette("Kriegeskorte_92_Images")` reproduces the published
human-IT model-RDM ranking on the 92-image stimulus set, and uses
`rdm_model_space_connectivity()` to show that within-subject sessions
cluster more tightly than between-subject pairs through the model space.
`vignette("Haxby_2001")` demonstrates the classification side
(`mvpa_model()` + `run_regional()` on the published VT mask reaches
91.7 % 8-way accuracy, recovering the canonical face/house pattern).

## 9. Relationship to prior work

The projection-onto-model-space view is implicit in the original RSA
framework of Kriegeskorte et al. (2008) and made explicit in
multivariate-RSA / multiple-RDM work (Khaligh-Razavi & Kriegeskorte
2014; Diedrichsen & Kriegeskorte 2017). What is novel in rMVPA is
treating the projection as a *first-class return value* of every RSA
fit — preserved through the regional and searchlight engines as a
list-column / attribute — and exposing the second-order operation
`F Fᵀ` as the canonical between-unit summary, with a memory-bounded
anchor extension for searchlight.

The framework relates to but is distinct from several adjacent ideas:

- **Informational connectivity** (Coutanche & Thompson-Schill 2013)
  asks whether decoding accuracies in two regions covary across trials.
  Model-space connectivity instead asks whether two regions share the
  *same* representational axes, irrespective of decoding accuracy.

- **Informational hierarchies via RSA** (Khaligh-Razavi et al. 2017;
  Cichy et al. 2016) use Spearman correlation of full RDM vectors
  across regions. That is recovered as the `raw_similarity` field in
  rMVPA's output (un-projected `Y Yᵀ / (p−1)`); the model-space
  projection is the additional, model-restricted view.

- **Feature-space prediction** (Mitchell et al. 2008; Anderson et al.
  2017) maps neural patterns to a feature space and quantifies
  decoding. `feature_rsa_connectivity()` extends this to ROI-to-ROI
  generalisation. Fingerprint-based connectivity is a complementary
  path that needs no fit and no feature space — only an explicit set of
  candidate model RDMs.

## 10. Limitations and open questions

The framework is honest about three commitments:

1. **Connectivity is restricted to the supplied model subspace.** If
   two regions share representational geometry that lies *outside* the
   model subspace, fingerprint connectivity won't see it. This is by
   design — it inherits the falsifiability of declared-model RSA — but
   readers should not interpret a low fingerprint similarity as
   evidence that two regions are unrelated. The `residual_similarity`
   field reports geometry shared *outside* the model subspace and
   serves as the appropriate complement.

2. **The basis depends on standardisation choice.** Pearson and
   Spearman variants of the basis differ; the choice should match the
   estimator used for the per-unit RSA scores. `model_space_connectivity()`
   exposes this via `method = "pearson" | "spearman"`.

3. **Anchor coverage is an empirical question.** k-means anchors with
   default `k = 20` are a reasonable starting point but do not
   guarantee that every meaningful brain mode is represented. Users
   should sweep `k` and inspect the eigenvalue spectrum of the cluster
   solution before drawing strong inferential conclusions from anchor
   maps.

Future work in the package will likely include (i) a leverage-score
anchor implementation with documented spectral approximation bounds,
(ii) an explicit JL-sketch storage path parallel to `return_fingerprint`,
and (iii) a permutation testing scaffold for inferential statements
on `F Fᵀ` entries.

## 11. Citation

If you use this framework in published work, please cite the package and
the foundational RSA literature:

> Buchsbaum, B. (2026). *rMVPA: Multivoxel Pattern Analysis in R*.
> R package version 0.1.2. <https://github.com/bbuchsbaum/rMVPA>

> Kriegeskorte N, Mur M, Bandettini P (2008). Representational
> similarity analysis — connecting the branches of systems neuroscience.
> *Frontiers in Systems Neuroscience* 2:4.

> Diedrichsen J, Kriegeskorte N (2017). Representational models: a
> common framework for understanding encoding, pattern-component, and
> representational-similarity analysis. *PLoS Computational Biology*
> 13(4): e1005508.

## References

Anderson AJ, Binder JR, Fernandino L, Humphries CJ, Conant LL, Aguilar
M, Wang X, Doko D, Raizada RDS (2017). Predicting neural activity
patterns associated with sentences using a neurobiologically motivated
model of semantic representation. *Cerebral Cortex* 27(9): 4379–4395.

Cichy RM, Khosla A, Pantazis D, Torralba A, Oliva A (2016). Comparison
of deep neural networks to spatio-temporal cortical dynamics of human
visual object recognition reveals hierarchical correspondence.
*Scientific Reports* 6:27755.

Coutanche MN, Thompson-Schill SL (2013). Informational connectivity:
identifying synchronized discriminability of multi-voxel patterns
across the brain. *Frontiers in Human Neuroscience* 7:15.

Diedrichsen J, Kriegeskorte N (2017). Representational models: a
common framework for understanding encoding, pattern-component, and
representational-similarity analysis. *PLoS Computational Biology*
13(4): e1005508.

Haxby JV, Gobbini MI, Furey ML, Ishai A, Schouten JL, Pietrini P (2001).
Distributed and overlapping representations of faces and objects in
ventral temporal cortex. *Science* 293: 2425–2430.

Khaligh-Razavi S-M, Kriegeskorte N (2014). Deep supervised, but not
unsupervised, models may explain IT cortical representation.
*PLoS Computational Biology* 10(11): e1003915.

Khaligh-Razavi S-M, Henriksson L, Kay K, Kriegeskorte N (2017).
Fixed versus mixed RSA: explaining visual representations by fixed and
mixed feature sets from shallow and deep computational models.
*Journal of Mathematical Psychology* 76: 184–197.

Kriegeskorte N, Mur M, Bandettini P (2008). Representational similarity
analysis — connecting the branches of systems neuroscience.
*Frontiers in Systems Neuroscience* 2:4.

Kriegeskorte N, Mur M, Ruff DA, Kiani R, Bodurka J, Esteky H, Tanaka K,
Bandettini PA (2008). Matching categorical object representations in
inferior temporal cortex of man and monkey. *Neuron* 60: 1126–1141.

Mitchell TM, Shinkareva SV, Carlson A, Chang K-M, Malave VL, Mason RA,
Just MA (2008). Predicting human brain activity associated with the
meanings of nouns. *Science* 320: 1191–1195.

Nili H, Wingfield C, Walther A, Su L, Marslen-Wilson W, Kriegeskorte N
(2014). A toolbox for representational similarity analysis. *PLoS
Computational Biology* 10(4): e1003553.
