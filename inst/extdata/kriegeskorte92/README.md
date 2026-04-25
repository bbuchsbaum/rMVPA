# Kriegeskorte 92-image dataset

`rdms.rds` is a compact bundle derived from the Kriegeskorte 2008 Neuron
data. It contains:

- **`brain`** — 8 human-IT representational dissimilarity matrices
  (4 subjects × 2 sessions, each 92×92, correlation-distance units).
- **`model`** — 8 model RDMs at the same 92-stimulus grid:
  *animacy*, *FaceBodyManmadeNatobj*, *monkeyIT*, *EVA*, *HMAX*, *V1*,
  *Silhouette*, *RADON*.
- **`source`** — citation of the original dataset.
- **`n_items`**, **`subjects`**, **`sessions`** — metadata.

## Source and license

Original publication:

> Kriegeskorte, Mur, Ruff, Kiani, Bodurka, Esteky, Tanaka, Bandettini (2008).
> *Matching categorical object representations in inferior temporal cortex of
> man and monkey.* **Neuron 60**:1126–1141.
> <https://doi.org/10.1016/j.neuron.2008.10.043>

The matrices in this bundle were re-derived from the supplemental MATLAB
files redistributed by the rsatoolbox project at
<https://github.com/rsagroup/rsatoolbox/tree/main/demos/92imageData>. We
bundle the two derived RDM files (`92_brainRDMs.mat`, `92_modelRDMs.mat`)
as a single XZ-compressed RDS to keep the package payload small.

When publishing analyses based on this bundle, cite Kriegeskorte et al.
(2008) and the rsatoolbox group.

## Regenerating the bundle

`data-raw/kriegeskorte92.R` is a one-shot R script that re-downloads the
upstream `.mat` files and rebuilds `rdms.rds` from scratch. Run it from the
package root if the upstream layout ever changes.
