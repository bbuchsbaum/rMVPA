# Select Features

Given a `feature_selection` specification object and a dataset, returns
the set of selected features as a binary vector.

## Usage

``` r
select_features(obj, X, Y, ...)

# S3 method for class 'catscore'
select_features(obj, X, Y, ranking.score = c("entropy", "avg", "max"), ...)

# S3 method for class 'FTest'
select_features(obj, X, Y, ...)
```

## Arguments

- obj:

  The `feature_selection` object specifying the feature selection method
  and its parameters.

- X:

  The dataset containing the training features. This can be a matrix or
  a `ROIVolume` or `ROISurface` object.

- Y:

  The dependent variable as a factor or numeric variable.

- ...:

  Additional arguments to be passed to the method-specific function.

- ranking.score:

  The feature score to use. Supported scores are "entropy", "avg", or
  "max". Default is "entropy".

## Value

A logical vector indicating the columns of `X` matrix that were
selected.

## Examples

``` r
fsel <- feature_selector("FTest", "top_k", 2)
coords <- rbind(c(1,1,1), c(2,2,2), c(3,3,3))
space <- neuroim2::NeuroSpace(c(10,10,10))
roi_data <- matrix(rnorm(100*3), 100, 3)
ROI <- neuroim2::ROIVec(space, coords=coords, roi_data)
Y <- factor(rep(c("a", "b"), each=50))
featureMask <- select_features(fsel, neuroim2::values(ROI), Y)
#> selecting features via FTest
#> cutoff type top_k
#> cutoff value 2
#> retaining 2 features in matrix with 3 columns
sum(featureMask) == 2
#> [1] TRUE

fsel2 <- feature_selector("FTest", "top_p", .1)
featureMask <- select_features(fsel2, neuroim2::values(ROI), Y)
#> selecting features via FTest
#> cutoff type top_p
#> cutoff value 0.1
#> retaining 1 features in matrix with 3 columns
```
