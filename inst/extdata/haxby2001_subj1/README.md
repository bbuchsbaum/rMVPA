# Haxby 2001 — subject 1 ventral-temporal patterns

`patterns.rds` is a small derivative of the Haxby et al. (2001) Subject 1
dataset, restricted to the published ventral-temporal (VT) mask.

## Contents

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

## Regenerating the bundle

`data-raw/haxby2001_subj1.R` is a one-shot R script that re-downloads the
upstream archive, applies the VT mask, computes per-(category, run) mean
patterns, and rebuilds `patterns.rds`. Run it from the package root if the
upstream layout changes.
