# rmvpa studio (web app)

React + Vite source for the rmvpadsl studio. `npm run build` writes the bundle
to `../rmvpadsl/inst/app/`, which the R package serves. `npm run build:example`
builds a static copy that reads the recorded example project in
`public/example/` (regenerate it with `Rscript fixture/make-fixture.R`; it
needs rmvpadsl installed). The example's project files, `explain` output and
`check` findings are real rmvpadsl output; its maps and numbers are synthetic.
