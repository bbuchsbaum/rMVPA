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
#> [1,] 0.4781820 0.5212485 5.694978e-04
#> [2,] 0.2373771 0.7626209 1.981426e-06
#> [3,] 0.3000222 0.6999670 1.087359e-05
#> [4,] 0.5335713 0.4646477 1.780954e-03
#> [5,] 0.1959692 0.8040303 5.324302e-07
#> [6,] 0.5248682 0.4736442 1.487679e-03
head(predict(fit, X, type = "scores"))
#>           [,1]
#> [1,] 0.7688035
#> [2,] 1.6677038
#> [3,] 1.4016334
#> [4,] 0.5820599
#> [5,] 1.8710969
#> [6,] 0.6116849
```
