# dsl/: the rmvpa analysis language

This folder is **not part of the rMVPA package**. It is listed in rMVPA's
`.Rbuildignore`, so CRAN builds of rMVPA never see it, and nothing in rMVPA
refers to it.

| Path | What it is |
|---|---|
| `rmvpadsl/` | R package: the YAML language, its compiler to rMVPA model specs, the run engine and the `rmvpa` command |
| `app/` | The studio's source (React + Vite). `npm run build` writes one bundle to `rmvpadsl/inst/app/`, which is committed, so installing rmvpadsl needs no Node. |

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

## The studio on a user's machine

The built app ships inside rmvpadsl (`inst/app/`). One bundle serves two ways:

* **`rmvpa studio`** / `rmvpadsl::rmvpa_studio()` serves it on 127.0.0.1 with
  httpuv and answers `/api/*` from the project's outputs. The printed address
  carries a session token; requests with another Host or Origin, or without the
  token, are refused.
* **`rmvpa export <dir>`** / `rmvpadsl::rmvpa_export_studio()` writes the app and
  every payload to a folder as `<script>` files. It opens from disk with no R
  and no server, or from any static host.

`.github/workflows/dsl.yaml` rebuilds the bundle and fails if the committed
`inst/app/` differs, then runs the rmvpadsl tests.
