# List the rMVPA API Lifecycle Registry

Returns a compact registry describing which exported entry points are
considered `"stable"`, `"experimental"`, or `"developer"`.

## Usage

``` r
rmvpa_api_lifecycle()
```

## Value

A data frame with columns `symbol`, `lifecycle`, and `notes`.
