# Proposal – Haar Pyramind Feature Extraction

The goal of this proposal is to incorporate a multi-scale Haar-wavelet pyramid ("Haar Pyramind") into the `rMVPA` package.  The transform should provide spatially localized features at multiple resolutions for use in searchlight and ROI analyses.

## 1. Background

Wavelet pyramids recursively decompose an image (or volumetric data) into low-pass and high-pass components.  A Haar basis is attractive due to its simplicity and the ability to exactly reconstruct the original signal from the pyramid coefficients.  Multi-resolution features can improve pattern analysis by capturing structure that spans different spatial scales.

## 2. Ticket HWT-S1-9 – Sprint Implementation

Ticket **HWT-S1-9** introduces normalization and coefficient bookkeeping for the Haar Pyramind.  Each level of the pyramid must retain enough metadata to accurately reconstruct voxel locations when inverting the transform.  This ensures compatibility with the existing searchlight infrastructure.

### Requirements

- Each decomposition level stores the scale factor, orientation (x, y, z), and voxel indices of the contributing source data.
- The reconstruction function must be able to combine coefficients from all levels to recreate the original pattern without loss.
- High-pass bands are zero-padded to the original voxel grid so that shapes remain consistent across levels.

### Proposed Implementation Steps

1. Add an S3 class `haar_pyramind` with fields:
   - `levels`: list of coefficient arrays (one per resolution level).
   - `meta`: metadata table containing scale factors and index mappings.
2. Implement `build_haar_pyramind(data, n_levels)` to recursively compute low/high-pass pairs.
3. Implement `reconstruct_from_pyramind(pyr)` to invert the transform using the stored metadata.
4. Integrate the pyramid generation with existing preprocessing so that `contrast_rsa_model` can accept pyramid coefficients as input features.

## 3. Open Questions

- Should optional smoothing be applied to the lowest-resolution level?
- How will noise whitening interact with the pyramid coefficients?

---

This document will evolve as development progresses.  The checklist below tracks completion of ticket HWT-S1-9 tasks.

- [ ] S3 class definition and constructor
- [ ] Forward decomposition (`build_haar_pyramind`)
- [ ] Reconstruction (`reconstruct_from_pyramind`)
- [ ] Integration with `contrast_rsa_model`
