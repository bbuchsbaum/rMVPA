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
#>                 a            b            c
#> [1,] 6.886596e-01 0.3113266184 1.374041e-05
#> [2,] 3.441882e-01 0.5970471490 5.876462e-02
#> [3,] 3.701903e-01 0.5928850519 3.692466e-02
#> [4,] 6.215954e-01 0.3783127781 9.183335e-05
#> [5,] 6.470987e-01 0.3528557085 4.562126e-05
#> [6,] 1.763322e-05 0.0002020792 9.997803e-01
head(predict(fit, X, type = "scores"))
#>            [,1]
#> [1,]  1.6436056
#> [2,]  0.0459737
#> [3,]  0.1408114
#> [4,]  1.2903446
#> [5,]  1.4208813
#> [6,] -2.1972109
```
