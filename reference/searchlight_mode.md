# Get Or Set Searchlight Mode

Simplified high-level runtime switch for searchlight performance
behavior. `"fast"` is the default optimized mode and `"legacy"` restores
conservative compatibility behavior.

## Usage

``` r
searchlight_mode(mode = NULL)
```

## Arguments

- mode:

  Optional mode value. One of `"fast"` or `"legacy"`. If `NULL`, returns
  current effective mode.

## Value

If `mode` is `NULL`, returns current mode as character. Otherwise, sets
the mode and returns it invisibly.
