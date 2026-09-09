# Display patterns in orthogonal spatial and target coordinates

Display patterns in orthogonal spatial and target coordinates

## Usage

``` r
rotate_patterns(fit, spatial = "varimax", target = "varimax")

# S3 method for class 'pattern_view'
predict(object, ...)
```

## Arguments

- fit:

  A fitted pattern model or global result with a refit.

- spatial, target:

  Orthogonal matrices, `"varimax"`, or `"none"`.

- object:

  A pattern view.

- ...:

  Arguments passed to `predict.pattern_fit`.

## Value

A `pattern_view` containing `L_b = A Q_b`, \\H = Q_b^T Q_t\\, and
`L_t = C Q_t`, in the fitted feature and whitened target coordinates.
Thus \\L_b H L_t^T\\ reconstructs \\A C^T\\. Back-transforms, a basis
identifier, and the Frobenius reconstruction error (evaluated without a
full feature by target matrix) are retained.

## Details

The underlying fit is unchanged. Predictions (including calibrated
scores) and invariant maps delegate to that fit and are bit-identical.
Display columns are not separately identified scientific components.
