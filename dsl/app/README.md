# rmvpa studio (web app)

React + Vite source for the rmvpadsl studio. `npm run build` writes the bundle
to `../rmvpadsl/inst/app/` (commit it: CI fails when it is stale). The same
bundle is used three ways, chosen at run time in `src/api.js`:

| Where the data comes from | When |
|---|---|
| `window.__RMVPA_DATA__`, set by `data/*.js` scripts | a static export (`rmvpa export <dir>`), opened from disk or a static host |
| `/api/*` with the `X-Rmvpa-Token` header | served by `rmvpa studio`, which puts the token in a `<meta>` tag |
| `./example/*.json` | `npm run dev` and the hosted example build |

The bundle is a single classic script (IIFE) loaded with `defer`, with CSS
injected by the script and the surface code inlined, because browsers refuse
ES-module scripts on `file://` pages. `npm run build:example`
builds a static copy that reads the recorded example project in
`public/example/` (regenerate it with `Rscript fixture/make-fixture.R`; it
needs rmvpadsl installed). The example's project files, `explain` output and
`check` findings are real rmvpadsl output; its maps and numbers are synthetic.

## Display libraries

- **surfview** (`ParcelSurface`-style rendering via `DataLayer` + `CurvatureLayer`)
  draws parcel values on the cortex (`src/surface.jsx`). The server offers a
  surface per regional analysis (exact for Schaefer-2018 from CBIG's
  annotation, approximate for other MNI atlases) and builds it when the user
  clicks Fetch (`src/assets.jsx`, `R/dsl_studio_assets.R`). The background
  picker works the same way for the image under the slices.
- **neuroatlasjs** (pinned to a commit; not on npm) holds atlas semantics:
  `VolumeAtlas` / `SurfaceAtlas` share one region table, and `ParcelValues`
  are projected to voxels and to vertices.
- Slices are painted on canvas (`src/brain.jsx`). neuroimjs's synchronized
  viewers are the next step where interactive orthogonal views are needed.
- neuromosaic (R) is the right home for exported, atlas-annotated static
  reports; it is not part of the live app.

`fixture/make-regions.R` builds the regional example from real files: CBIG
Schaefer-2018 (200 parcels, 7 networks; MNI volume and fsaverage5
annotations) and TemplateFlow fsaverage5 pial + sulc. Accuracies are
synthetic.
