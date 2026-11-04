# dsl/: the rmvpa analysis language

This folder is **not part of the rMVPA package**. It is listed in rMVPA's
`.Rbuildignore`, so CRAN builds of rMVPA never see it, and nothing in rMVPA
refers to it.

| Path | What it is |
|---|---|
| `rmvpadsl/` | R package: the YAML language, its compiler to rMVPA model specs, the run engine and the `rmvpa` command |
| `app/` | Reserved for the web UI's source (React + Vite). Its build output goes to `rmvpadsl/inst/app/`, served by the R package. |

## Boundary

* rmvpadsl uses only rMVPA's **exported** API (`@import rMVPA`). It never
  calls `rMVPA:::` internals.
* Facts about rMVPA's behaviour that the language encodes (which scopes each
  model supports, which permutation nulls exist, which optional packages a
  verb needs) live in `.dsl_scope_support()`, `.dsl_inference_support()` and
  `.dsl_required_packages()`. The end-to-end tests run every verb against the
  real constructors, so a change in rMVPA that breaks one of these facts fails
  the rmvpadsl tests.

## Developing

```r
# from the repository root
install.packages(".", repos = NULL, type = "source")    # rMVPA
devtools::test("dsl/rmvpadsl")
```

The command can be run from a checkout with `dsl/rmvpadsl/exec/rmvpa help`
(needs pkgload), or installed with `rmvpadsl::install_rmvpa_command(dir)`.
