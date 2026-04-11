# Merge Multiple Results

Merge Multiple Results

## Usage

``` r
merge_results(obj, result_set, indices, id, ...)
```

## Arguments

- obj:

  The base object containing merge specifications

- result_set:

  List of results to be merged

- indices:

  List of indices corresponding to each result

- id:

  Identifier for the merged result

- ...:

  Additional arguments passed to specific merge methods

## Value

A merged result object containing:

- Combined results from all input objects

- Associated indices

- Merged metadata

## Examples

``` r
# \donttest{
ds <- gen_sample_dataset(D = c(5, 5, 5), nobs = 20, nlevels = 2)
model <- load_model("sda_notune")
mspec <- mvpa_model(
  model   = model,
  dataset = ds$dataset,
  design  = ds$design,
  model_type = "classification"
)

# Construct a minimal result_set resembling a single CV fold:
obs  <- ds$design$y_train
levs <- levels(obs)
prob_mat <- matrix(
  1 / length(levs),
  nrow = length(obs),
  ncol = length(levs),
  dimnames = list(NULL, levs)
)

result_set <- tibble::tibble(
  test_ind      = list(seq_along(obs)),
  probs         = list(prob_mat),
  error         = FALSE,
  error_message = "~"
)

merge_results(mspec, result_set, indices = list(1:10), id = 1)
#> # A tibble: 1 × 6
#>   result         indices    performance    id error error_message
#>   <list>         <list>     <list>      <dbl> <lgl> <chr>        
#> 1 <bnry_cl_ [6]> <list [1]> <dbl [2]>       1 FALSE ~            
# }
```
