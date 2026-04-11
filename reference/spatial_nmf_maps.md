# Spatial NMF on Map Lists

**Raison d'etre.** Spatial NMF factorizes a non-negative
subject-by-voxel matrix \\X\\ into \\W\\ (subject loadings) and \\H\\
(spatial components) such that \\X \approx W H\\. \\W\\ captures how
strongly each subject expresses each component, while \\H\\ encodes
spatially interpretable, additive patterns. An optional graph Laplacian
penalty encourages smooth, neuroanatomically plausible components. The
resulting representation is compact and supports downstream inference
(component tests, global CV tests, and stability).

## Usage

``` r
spatial_nmf_maps(
  group_A,
  group_B = NULL,
  mask = NULL,
  dims = NULL,
  k,
  lambda = 0,
  fast = FALSE,
  graph = NULL,
  neighbors = 6,
  na_action = c("zero", "error"),
  return_maps = TRUE,
  return_data = FALSE,
  component_test = NULL,
  global_test = NULL,
  stability = NULL,
  voxelwise_stats = c("none", "stability_zp"),
  parallel = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- group_A:

  List of NeuroVol/NeuroSurface maps for group A. All maps must share
  the same spatial grid/geometry; if maps carry sparse indices, those
  indices must match.

- group_B:

  Optional list of NeuroVol/NeuroSurface maps for group B (same
  requirements as group_A).

- mask:

  Mask object for volumetric data (NeuroVol or logical/numeric vector).
  Required for volumetric inputs; optional for surface inputs.

- dims:

  Optional spatial dimensions for volumetric masks given as vectors.

- k:

  Number of components.

- lambda:

  Spatial regularization strength (0 = none).

- fast:

  Logical; use faster, lower-iteration defaults in the NMF fit (e.g.,
  random init, fewer iterations). You can still override specific
  optimization settings via \`...\`.

- graph:

  Optional graph Laplacian list (from build_graph_laplacian). If
  \`graph\$A\` is weighted, set \`graph\$weighted=TRUE\` to preserve
  weights; otherwise edges are binarized.

- neighbors:

  Neighborhood size for volumetric adjacency (6/18/26).

- na_action:

  How to handle NA values in maps: "zero" or "error".

- return_maps:

  Logical; return component maps as NeuroVol/NeuroSurface.

- return_data:

  Logical; include the data matrix in the result.

- component_test:

  NULL to skip, TRUE for defaults, or a list of arguments passed to
  spatial_nmf_component_test.

- global_test:

  NULL to skip, TRUE for defaults, or a list of arguments passed to
  spatial_nmf_global_test.

- stability:

  NULL to skip, TRUE for defaults, or a list of arguments passed to
  spatial_nmf_stability.

- voxelwise_stats:

  Optional voxelwise statistics to compute. Use "stability_zp" to return
  bootstrap z- and p-value maps derived from \`stability\` (forces
  \`stability\$return_maps = TRUE\`).

- parallel:

  Logical; enable parallel processing for inference functions
  (component_test, global_test, stability). If NULL (default),
  auto-detects based on active future plan and availability of
  future.apply.

- progress:

  Logical; report progress via progressr package (works with parallel
  futures). Default FALSE.

- ...:

  Additional arguments passed to spatial_nmf_fit.

## Value

A list with fields:

- fit: spatial_nmf_fit object.

- components: list of spatial component maps (if return_maps).

- groups: factor of group labels for each subject.

- mask_indices: indices used to vectorize maps.

- map_type: "volume" or "surface".

- data: subject-by-voxel matrix (if return_data).

- voxelwise: list of voxelwise z/p maps (if requested).

## Details

**Why not just voxelwise t-tests?** Spatial NMF targets multivariate
structure by learning *patterns of covarying voxels* rather than testing
each voxel independently. This yields a low-dimensional representation
that is often more interpretable, statistically powerful (fewer multiple
comparisons), and robust to noise. With optional spatial regularization,
components are anatomically smoother than voxelwise maps, and the
framework provides component-level inference and stability analyses that
reflect network-level effects.

**Interpreting inference outputs.**

- *Component tests* (\`spatial_nmf_component_test\`): use \`p_fwer\` to
  identify components whose loadings differ reliably between groups;
  \`p_unc\` is uncorrected. \`p_global\` tests whether any component
  differs.

- *Global test* (\`spatial_nmf_global_test\`): \`stat\` is
  cross-validated performance (AUC or accuracy) on component loadings; a
  significant \`p_value\` indicates a multivariate group difference
  without pinpointing a component.

- *Stability* (\`spatial_nmf_stability\`): interpret \`selection\` as
  voxelwise stability (frequency among top-loadings), and \`cv\` as
  variability; higher \`component_similarity\` indicates stable
  component identity.

Convenience wrapper that converts lists of NeuroVol/NeuroSurface maps
into a subject-by-voxel matrix and fits spatially regularized NMF.

If \`component_test\`, \`global_test\`, or \`stability\` are requested
and their argument lists do not specify \`parallel\`,
\`spatial_nmf_maps\` will enable parallel execution automatically when a
future plan with more than one worker is active and the \`future.apply\`
package is available. If \`fast = TRUE\`, the same setting is passed
through to \`global_test\` and \`stability\` (unless explicitly
overridden) so bootstrap/CV refits use the faster defaults too.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Requires list of NeuroVol objects per group
  # result <- spatial_nmf_maps(group_A_vols, group_B_vols, k=5)
} # }
```
