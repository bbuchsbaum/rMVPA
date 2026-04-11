# Merge Predictions

Combine predictions from multiple models on the same test set.

## Usage

``` r
merge_predictions(obj1, rest, ...)

# S3 method for class 'regression_prediction'
merge_predictions(obj1, rest, ...)

# S3 method for class 'classification_prediction'
merge_predictions(obj1, rest, ...)
```

## Arguments

- obj1:

  The first object containing predictions.

- rest:

  Other objects containing predictions.

- ...:

  Additional arguments. Methods for this generic may implement specific
  arguments such as \`weights\` to control how predictions are combined.

## Value

A combined object with merged predictions.

## Examples

``` r
# \donttest{
p1 <- structure(
  list(class = c("a","b"),
       probs = matrix(c(.8,.2,.3,.7), ncol=2, dimnames=list(NULL,c("a","b")))),
  class = c("classification_prediction", "prediction", "list")
)
p2 <- structure(
  list(class = c("b","a"),
       probs = matrix(c(.4,.6,.6,.4), ncol=2, dimnames=list(NULL,c("a","b")))),
  class = c("classification_prediction", "prediction", "list")
)
merge_predictions(p1, list(p2))
#> $class
#> [1] "a" "b"
#> 
#> $probs
#>        a    b
#> [1,] 0.6 0.45
#> [2,] 0.4 0.55
#> 
#> attr(,"class")
#> [1] "classification_prediction" "prediction"               
#> [3] "list"                     
# }
```
