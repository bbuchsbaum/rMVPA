# Clean Up Shared-Memory Segments

Closes shared-memory segments created by
[`use_shard`](http://bbuchsbaum.github.io/rMVPA/reference/use_shard.md).
Called automatically at the end of
[`mvpa_iterate()`](http://bbuchsbaum.github.io/rMVPA/reference/mvpa_iterate.md);
can also be invoked manually.

## Usage

``` r
shard_cleanup(shard_data)
```

## Arguments

- shard_data:

  List returned by `shard_prepare_dataset` (i.e.\\
  `mod_spec$shard_data`).

## Value

`NULL` (invisibly).
