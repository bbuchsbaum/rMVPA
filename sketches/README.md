# sketches/

Exploratory, **non-shipping** material. This whole directory is excluded from
the package build via `.Rbuildignore`.

## rsa_schematic.R

A prototype generator for two RSA pipeline schematics:

- `rsa_schematic_basic()` -- standard RSA. Five panels:
  1. Patterns (conditions × voxels) heatmap
  2. Neural RDM (with `dist()` arrow)
  3a/3b. Two model RDMs (e.g. category, identity)
  4. Lower-triangle scatter of model vs. neural with Spearman ρ
  5. Multi-RDM regression β bar plot

- `rsa_schematic_msreve()` -- contrast / MS-ReVE RSA. Six panels:
  1. Contrast matrix `C` (K × Q)
  2a/2b. Predicted Δ_q = c_q c_qᵀ for each contrast
  3. Neural RDM
  4. β_q bar plot
  5. Signed voxel map β_q · Δ_{q,v}
  6. Equation/legend

Run it:

```r
source("sketches/rsa_schematic.R")
rsa_schematic_basic("sketches/rsa_basic.png")
rsa_schematic_msreve("sketches/rsa_msreve.png")
```

Both functions also return the underlying matrices invisibly so they can be
re-used in a vignette chunk without re-computing.

## Status

Prototype. Not wired into any vignette yet -- still deciding whether a static
schematic is more useful than the current per-step exposition.
