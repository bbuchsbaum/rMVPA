# Schema-driven combiner for searchlight results

Drop-in replacement for
[`combine_standard`](http://bbuchsbaum.github.io/rMVPA/reference/combine_standard.md)
that uses the model's
[`output_schema`](http://bbuchsbaum.github.io/rMVPA/reference/output_schema.md)
to build output maps. Falls back to `combine_standard` if the schema is
`NULL`.

## Usage

``` r
combine_schema_standard(model_spec, good_results, bad_results)
```

## Arguments

- model_spec:

  A list containing the model specification

- good_results:

  A data frame containing the successful classifier results

- bad_results:

  A data frame containing the unsuccessful classifier results

## Value

A `searchlight_result` object.
