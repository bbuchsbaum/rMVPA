# Compute Group Means of a Matrix

This function calculates the average vector for each level of a grouping
variable in a given matrix.

## Usage

``` r
group_means(X, margin, group)
```

## Arguments

- X:

  A matrix for which group means should be calculated.

- margin:

  An integer specifying the margin to average over. Use 1 for averaging
  over rows, and 2 for averaging over columns.

- group:

  A grouping variable, either a factor or an integer vector, that
  defines the groups to calculate the means for.

## Value

A matrix with the same number of rows or columns (depending on the
margin) as the input matrix X, and the number of columns or rows
corresponding to the number of unique groups in the grouping variable.

## Examples

``` r
# Create a random matrix
data <- matrix(rnorm(100 * 100), 100, 100)

# Define a grouping variable
groups <- factor(rep(1:5, each = 20))

# Calculate group means for each row
row_means <- group_means(data, margin = 1, group = groups)

# Calculate group means for each column
col_means <- group_means(data, margin = 2, group = groups)
```
