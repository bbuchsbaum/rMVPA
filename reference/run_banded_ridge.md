# Run a chunked whole-brain banded-ridge encoding analysis

\`run_banded_ridge()\` evaluates each active mask voxel exactly once.
Optional \`response_partitions\` may reorder non-overlapping response
batches (the same contract used by regional batching), but may not
overlap or omit voxels. When \`model\$delta_sets\` is non-empty, the
result additionally contains \`predictive_leave_one_band_out\`: matched
full/reduced OOF metrics and predictions, independently selected reduced
hyperparameters, spatial \`delta_cv_r2\_\<band\>\` maps, and B+1
model-cost provenance.

## Usage

``` r
run_banded_ridge(
  model,
  target_batch_size = model$target_batch_size,
  response_partitions = NULL
)
```

## Arguments

- model:

  A \`banded_ridge_model\` specification.

- target_batch_size:

  Optional runtime override of the maximum response chunk size.

- response_partitions:

  Optional list of non-overlapping global voxel IDs that exactly cover
  the active mask.

## Value

A \`banded_ridge_result\` with spatial maps, metrics, exact outer-fold
hyperparameters, optional predictions/weights, and allocation
provenance.
