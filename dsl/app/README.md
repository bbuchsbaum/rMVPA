# rmvpa studio (web app)

React + Vite source for the rmvpadsl studio. `npm run build` writes the bundle
to `../rmvpadsl/inst/app/`, which the R package serves. `npm run build:example`
builds a static copy that reads the recorded example project in
`public/example/` (regenerate it with `Rscript fixture/make-fixture.R`; it
needs rmvpadsl installed). The example's project files, `explain` output and
`check` findings are real rmvpadsl output; its maps and numbers are synthetic.

## Display libraries

- **surfview** (`ParcelSurface`-style rendering via `DataLayer` + `CurvatureLayer`)
  draws parcel values on the cortex. It is loaded lazily (`src/surface.jsx`),
  so three.js stays out of the main bundle.
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
