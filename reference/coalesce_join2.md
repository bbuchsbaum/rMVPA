# Coalesce Join Two Data Frames

This function performs a specified type of join on two data frames and
then coalesces the joined columns based on their common column names.

## Usage

``` r
coalesce_join2(
  x,
  y,
  by = NULL,
  suffix = c(".x", ".y"),
  join = dplyr::full_join,
  ...
)
```

## Arguments

- x:

  A data frame to be joined.

- y:

  A second data frame to be joined.

- by:

  A character vector of variables to join by. If NULL (the default), the
  function will use the common column names in 'x' and 'y'.

- suffix:

  A character vector of length 2, specifying the suffixes to be used for
  making unique the common column names in 'x' and 'y'. The default is
  c(".x", ".y").

- join:

  A join function to be used for joining the data frames. The default is
  dplyr::full_join.

- ...:

  Additional arguments passed on to the join function.

## Value

A data frame resulting from the specified join operation and coalesced
columns.
