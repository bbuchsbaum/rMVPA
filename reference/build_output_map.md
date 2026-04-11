# Build Spatial Output Map for a Single Metric

Constructs a spatial object (e.g., `NeuroVol`, `NeuroSurface`) from a
numeric vector of per-center metric values and the corresponding center
IDs.

## Usage

``` r
build_output_map(dataset, metric_vector, ids, ...)
```

## Arguments

- dataset:

  The dataset object.

- metric_vector:

  Numeric vector of metric values (one per center).

- ids:

  Integer vector of center IDs corresponding to `metric_vector`.

- ...:

  Additional arguments.

## Value

A spatial object (`NeuroVol`, `NeuroSurface`, etc.)
