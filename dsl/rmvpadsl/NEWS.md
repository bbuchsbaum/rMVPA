# rmvpadsl (development version)

* The studio ships in the package (`inst/app`, one bundle built from
  `dsl/app`). `rmvpa_studio()` (`rmvpa studio`) serves it on 127.0.0.1 with a
  session token and Host/Origin checks; `rmvpa_export_studio()`
  (`rmvpa export <dir>`) writes it with the project's results to a folder
  that opens from disk with no R or server. Both read rmvpa_run() outputs only.
* Searchlight permutation runs also write the pooled null as 1001 quantiles
  (`*_stat-<metric>_null.tsv`), so its distribution can be redrawn.

* First version of the rmvpa analysis language (DSL 0.1), moved out of rMVPA
  into its own package so that rMVPA's CRAN package stays free of it.
  Verbs: decode, transfer, rsa, rsa_vector, rsa_contrast, rsa_feature,
  rsa_era, rsa_connect, rsa_mediate, encode, pattern, manova and custom.
  Commands: new, add, explain, check, draft, run, status, lock and migrate.
