# Wrap output results

This function wraps the output results of the performance matrix into a
list of spatial objects (NeuroVol or NeuroSurface) for each column in
the performance matrix, and structures it as a searchlight_result.

## Usage

``` r
wrap_out(perf_mat, dataset, ids = NULL)
```

## Arguments

- perf_mat:

  A performance matrix (voxels/vertices x metrics) containing classifier
  results.

- dataset:

  A dataset object containing the dataset information (including mask
  and type).

- ids:

  An integer vector of voxel/vertex indices corresponding to the rows of
  \`perf_mat\`. These are typically global indices into the mask space
  for volumetric data, or vertex numbers for surface data.

## Value

A \`searchlight_result\` object containing

- \`results\`: Named spatial maps for each metric.

- \`n_voxels\`: Total number of voxels/vertices defined by the mask.

- \`active_voxels\`: Number of voxels/vertices with results.

- \`metrics\`: Character vector of metric names.
