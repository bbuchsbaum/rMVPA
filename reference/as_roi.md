# Convert object to ROI

Convert the provided object into an ROIVolume or ROISurface object.

## Usage

``` r
as_roi(obj, data, ...)
```

## Arguments

- obj:

  The object to be converted.

- data:

  The associated data object.

- ...:

  Additional arguments passed to methods.

## Value

An ROIVolume or ROISurface object.

## Examples

``` r
# \donttest{
  ds <- gen_sample_dataset(c(5,5,5), 20, blocks=2)
  vox <- sample(which(ds$dataset$mask > 0), 10)
  samp <- data_sample(ds$dataset, vox)
  roi <- as_roi(samp, ds$dataset)
# }
```
