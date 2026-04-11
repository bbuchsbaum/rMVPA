# Create a feature selection specification

This function creates a feature selection specification using the
provided method, cutoff type, and cutoff value.

## Usage

``` r
feature_selector(method, cutoff_type, cutoff_value)
```

## Arguments

- method:

  The type of feature selection method to use. Supported methods are
  "FTest" and "catscore".

- cutoff_type:

  The type of threshold used to select features. Supported cutoff types
  are "top_k" and "top_p".

- cutoff_value:

  The numeric value of the threshold cutoff.

## Value

A list with a class name equal to the `method` argument.

## Details

The available feature selection methods are: - FTest: Computes a one-way
ANOVA for every column in the feature matrix. - catscore: Computes a
correlation adjusted t-test for every column in the matrix using
`sda.ranking` from the `sda` package.

## Examples

``` r
fsel <- feature_selector("FTest", "top_k", 1000)
fsel <- feature_selector("FTest", "top_p", .1)
class(fsel) == "FTest"
#> [1]  TRUE FALSE FALSE
```
