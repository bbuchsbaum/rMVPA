# Create a Surface-Based MVPA Dataset Object

Creates a dataset object for surface-based MVPA analysis that
encapsulates a training dataset, an optional test dataset, and a vertex
mask.

## Usage

``` r
mvpa_surface_dataset(train_data, test_data = NULL, mask = NULL, name = "")
```

## Arguments

- train_data:

  The training data set: must inherit from `NeuroSurfaceVector`

- test_data:

  Optional test data set: must inherit from `NeuroSurfaceVector`
  (default: NULL)

- mask:

  Optional binary mask for vertices. If NULL, creates mask from training
  data indices

- name:

  Optional label to identify the dataset (e.g., "lh" or "rh" to indicate
  hemisphere)

## Value

An `mvpa_surface_dataset` object (S3 class) containing:

- train_data:

  The training data as a `NeuroSurfaceVector` instance

- test_data:

  The test data as a `NeuroSurfaceVector` instance (if provided)

- mask:

  A numeric vector indicating valid vertices (1) and excluded vertices
  (0)

- name:

  Character string identifier for the dataset

- has_test_set:

  Logical flag indicating whether this dataset has a test set

## Details

If no mask is provided, one will be created automatically using the
indices from the training data. The mask will be a numeric vector with
length equal to the number of nodes in the surface geometry.

## See also

[`mvpa_dataset`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_dataset.md)
for creating volume-based MVPA datasets

[`mvpa_design`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_design.md)
for creating the corresponding design object

## Examples

``` r
if (FALSE) { # \dontrun{
# Create surface dataset with automatic mask
train_surf <- NeuroSurfaceVector(geometry, data)
dataset <- mvpa_surface_dataset(train_surf, name="lh")

# Create dataset with test data and custom mask
test_surf <- NeuroSurfaceVector(geometry, test_data)
mask <- numeric(length(nodes(geometry)))
mask[roi_indices] <- 1
dataset <- mvpa_surface_dataset(train_surf, test_surf, mask, name="rh")
} # }
```
