# Register a Custom MVPA Model

Adds a user-defined model specification to the rMVPA model registry
(MVPAModels).

## Usage

``` r
register_mvpa_model(name, model_spec)
```

## Arguments

- name:

  A character string, the unique name for the model.

- model_spec:

  A list containing the model specification. It must include elements:
  \`type\` ("Classification" or "Regression"), \`library\` (character
  vector of required packages for the \*model itself\*, not for rMVPA's
  wrappers), \`label\` (character, usually same as name), \`parameters\`
  (data.frame of tunable parameters: parameter, class, label), \`grid\`
  (function to generate tuning grid, takes x, y, len args), \`fit\`
  (function), \`predict\` (function), and \`prob\` (function for
  classification, takes modelFit, newdata; should return matrix/df with
  colnames as levels).

## Value

Invisibly returns the registered model specification.

## Examples

``` r
# A minimal spec: a nearest-centroid classifier using only base R
my_centroid_spec <- list(
  type = "Classification", library = NULL, label = "my_centroid",
  parameters = data.frame(parameter = "parameter", class = "character",
                          label = "parameter"),
  grid = function(x, y, len = NULL) data.frame(parameter = "none"),
  fit = function(x, y, wts, param, lev, last, weights, classProbs, ...) {
    list(centroids = rowsum(as.matrix(x), y) / as.vector(table(y)),
         levels = levels(y))
  },
  predict = function(modelFit, newdata, ...) {
    p <- cor(t(as.matrix(newdata)), t(modelFit$centroids))
    factor(modelFit$levels[max.col(p)], levels = modelFit$levels)
  },
  prob = function(modelFit, newdata, ...) {
    p <- exp(cor(t(as.matrix(newdata)), t(modelFit$centroids)))
    p <- p / rowSums(p)
    colnames(p) <- modelFit$levels
    p
  }
)
register_mvpa_model("my_centroid", my_centroid_spec)
mod <- load_model("my_centroid")
mod$label
#> [1] "my_centroid"
```
