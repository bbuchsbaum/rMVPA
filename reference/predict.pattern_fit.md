# Predict from a fitted pattern model

Predict from a fitted pattern model

## Usage

``` r
# S3 method for class 'pattern_fit'
predict(
  object,
  newdata,
  type = c("prob", "class", "decode", "scores", "encode"),
  targets = NULL,
  ...
)
```

## Arguments

- object:

  A `pattern_fit`.

- newdata:

  Numeric matrix of observations (rows) by features (columns): either
  all input features or only the retained ones.

- type:

  What to return: `"prob"` (class posterior probabilities; categorical
  targets), `"class"` (the most probable class), `"decode"`
  (posterior-mean targets on the original scale; for a categorical fit
  these are posterior-mean one-hot codes, which are not probabilities
  and do not sum to one, so use `"prob"` instead), `"scores"`
  (calibrated component scores \\z = G^{+} u\\), or `"encode"`
  (predicted brain measurements for supplied `targets`; returned over
  the retained features, with the retained column positions in attribute
  `"feature_index"`).

- targets:

  Targets for `type = "encode"`: a factor or a numeric vector/matrix on
  the original scale.

- ...:

  Ignored.

## Value

A matrix (or factor for `type = "class"`).

## Details

Classification uses the Gaussian class-conditional model implied by the
fit: \\p(c \mid x) \propto \pi_c \exp(m_c' u - m_c' G m_c / 2)\\ with
\\m_c = C' y_w(c)\\. Decoding uses the working prior \\y_w \sim N(0,
I)\\ on the whitened targets, giving \\\hat t = \Phi (I + G \Phi)^{-1}
u\\ and \\\hat y_w = C \hat t\\. Neither forms \\G^{-1}\\; a fit with no
retained signal returns the class priors or the target means.

## Examples

``` r
ds <- gen_sample_dataset(c(6, 6, 4), 60, nlevels = 3, blocks = 3)
spec <- pattern_model(ds$dataset, ds$design, rank = 1, refit = TRUE)
fit <- run_global(spec, refit = TRUE)$refit
X <- get_feature_matrix(ds$dataset)
head(predict(fit, X, type = "prob"))
#>              a         b            c
#> [1,] 0.4888936 0.5104231 6.832907e-04
#> [2,] 0.2537751 0.7462228 2.110850e-06
#> [3,] 0.3142492 0.6857393 1.143648e-05
#> [4,] 0.5400232 0.4578627 2.114068e-03
#> [5,] 0.2054227 0.7945769 4.399894e-07
#> [6,] 0.5171007 0.4816278 1.271492e-03
head(predict(fit, X, type = "scores"))
#>           [,1]
#> [1,] 0.7366898
#> [2,] 1.6581381
#> [3,] 1.3927199
#> [4,] 0.5514731
#> [5,] 1.9021100
#> [6,] 0.6351006
```
