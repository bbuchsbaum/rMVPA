# Build spatial output maps from an output schema

Uses a named schema list to select columns from a performance matrix and
create one spatial map per metric.

## Usage

``` r
build_maps_from_schema(schema, perf_mat, dataset, ids)
```

## Arguments

- schema:

  Named list; names are metric labels, values are `"scalar"` or
  `"vector[N]"`.

- perf_mat:

  Numeric matrix with rows = ROIs and columns = metrics.

- dataset:

  The dataset object (used by
  [`build_output_map`](http://bbuchsbaum.github.io/rMVPA/reference/build_output_map.md)).

- ids:

  Integer vector of ROI center IDs corresponding to rows of `perf_mat`.

## Value

A named list of spatial map objects.
