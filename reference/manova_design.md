# Create a MANOVA Design

This function creates a MANOVA design object containing a formula
expression and a named list of data.

## Usage

``` r
manova_design(formula, data)
```

## Arguments

- formula:

  A formula expression specifying the MANOVA regression model.

- data:

  A named list containing the dissimilarity matrices and any other
  auxiliary variables.

## Value

A MANOVA design object with class attributes "manova_design" and "list".

## Details

The function takes a formula expression and a named list of data as
input, and returns a MANOVA design object. The object is a list that
contains the formula expression and the named list of data with class
attributes "manova_design" and "list". This object can be further used
for MANOVA analysis or other related multivariate statistical methods.

## Examples

``` r
# Create simple dissimilarity matrices
dissimilarity_matrix_y  <- matrix(rnorm(9), nrow = 3)
dissimilarity_matrix_x1 <- matrix(rnorm(9), nrow = 3)
dissimilarity_matrix_x2 <- matrix(rnorm(9), nrow = 3)

# Create a MANOVA design
formula   <- y ~ x1 + x2
data_list <- list(
  y  = dissimilarity_matrix_y,
  x1 = dissimilarity_matrix_x1,
  x2 = dissimilarity_matrix_x2
)
manova_design_obj <- manova_design(formula, data_list)
```
