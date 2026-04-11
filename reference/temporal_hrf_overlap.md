# HRF-overlap temporal confound RDM

Builds a temporal confound RDM based on predicted HRF overlap between
events. Each event is convolved with a canonical HRF and pairwise
similarity is computed as either normalized overlap (cosine similarity)
or Pearson correlation. The result can be returned as a distance-like
quantity suitable for regression.

## Usage

``` r
temporal_hrf_overlap(
  onsets,
  durations = NULL,
  run = NULL,
  TR,
  hrf = c("spm", "glover", "gamma"),
  oversampling = 16L,
  length_s = 32,
  similarity = c("overlap", "corr"),
  metric = c("distance", "similarity"),
  within_blocks_only = TRUE,
  normalize = c("z", "rank", "none"),
  as_dist = TRUE
)
```

## Arguments

- onsets:

  numeric vector of event onsets in seconds

- durations:

  optional numeric vector of event durations in seconds (scalar or
  length N). Default 0 (impulse).

- run:

  optional vector of run/block identifiers

- TR:

  repetition time in seconds (required)

- hrf:

  character: one of `"spm"`, `"glover"`, or `"gamma"`

- oversampling:

  integer oversampling factor relative to TR (default 16)

- length_s:

  numeric length of HRF kernel in seconds (default 32)

- similarity:

  one of `"overlap"` (cosine similarity) or `"corr"`

- metric:

  return type: `"distance"` or `"similarity"` (default `"distance"`)

- within_blocks_only:

  logical; if TRUE zero-out cross-run entries (default TRUE)

- normalize:

  one of `"rank"`, `"z"`, or `"none"`

- as_dist:

  logical; if TRUE return a `dist` object (default TRUE)

## Value

A `dist` object or symmetric matrix (N x N)

## Examples

``` r
if (FALSE) { # \dontrun{
  onsets <- c(0, 5, 10, 15, 20)
  hrf_rdm <- temporal_hrf_overlap(onsets, TR=2, hrf="spm")
} # }
```
