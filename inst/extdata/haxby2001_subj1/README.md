# Haxby 2001 — subject 1 ventral-temporal patterns

This directory ships two derivatives of the Haxby et al. (2001) Subject 1
dataset, both restricted to the published ventral-temporal (VT) mask:

- **`patterns.rds`** (≈100 KB) — per-(category, run) **block-mean** patterns.
  Use this for the standard `mvpa_model()` / `rsa_model()` workflow.
- **`trials.rds`** (≈462 KB) — **trial-level** per-TR patterns. Use this
  when you need within-block residuals (the prerequisite for multivariate
  noise normalisation / MVNN), or for any analysis that benefits from the
  unaveraged trial-level structure.

Both are derived from the same upstream archive (PyMVPA's curated Subject 1
tutorial bundle) and carry identical voxel ordering and mask metadata.

## Contents of `patterns.rds`

- **`patterns`** — `[n_obs × n_voxels]` numeric matrix of per-(category, run)
  mean BOLD signals over the VT voxels. `n_obs = 8 × 12 = 96` (8 non-rest
  categories × 12 runs). `n_voxels = 577` (size of the VT mask).
- **`category`** — `factor` of length 96 with levels
  `bottle / cat / chair / face / house / scissors / scrambledpix / shoe`.
- **`run`** — `integer` of length 96, the original run/chunk index (0–11).
  Use this as the `block_var` for blocked cross-validation.
- **`mask_dim`**, **`mask_idx`**, **`mask_space`** — metadata describing how
  to back-project the 577 voxel columns into the original 40 × 64 × 64
  scanner space.
- **`source`** — citation of the original dataset.

## Contents of `trials.rds`

Same structure as `patterns.rds`, but each row is one *task* TR rather than
a block-mean:

- **`patterns`** — `[864 × 577]` numeric matrix.
- **`category`** — `factor` of length 864.
- **`run`** — `integer` of length 864 (12 unique values).
- **`block`** — `integer` of length 864. A "block" is a maximal run of
  consecutive same-category TRs within one run; there are 96 blocks total
  (8 per run × 12 runs). Use this as the unit for MVNN residual-covariance
  estimation: per block, subtract the block's mean pattern from each TR to
  obtain trial-level noise residuals.
- Same `mask_dim` / `mask_idx` / `mask_space` / `source` metadata as above.

## Source and license

Original publication:

> Haxby, Gobbini, Furey, Ishai, Schouten, Pietrini (2001).
> *Distributed and overlapping representations of faces and objects in
> ventral temporal cortex.* **Science 293**: 2425–2430.

The matrix in this bundle was derived from the PyMVPA-curated Subject 1
tutorial archive at
<http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz>
(`bold.nii.gz`, `mask4_vt.nii.gz`, `labels.txt`). We bundle only the
mean-pattern matrix (≈100 KB) rather than the raw 4D BOLD (≈300 MB) so the
package payload stays small. The original data is freely available for
non-commercial research; cite Haxby et al. (2001) when publishing.

## Regenerating the bundles

Two one-shot R scripts in `data-raw/`:

- `haxby2001_subj1.R` rebuilds `patterns.rds` (block-mean patterns) from the
  upstream archive.
- `haxby2001_subj1_trials.R` rebuilds `trials.rds` (per-TR patterns + block
  indices) from the same archive.

Both run from the package root and write directly under `inst/extdata/`.
