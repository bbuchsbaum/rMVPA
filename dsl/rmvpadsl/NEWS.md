# rmvpadsl (development version)

* Regional analyses get permutation tests when rMVPA has
  `run_permutation_regional()`: decode and RSA over regions default to
  `inference: { test: permutation }`, as the searchlight does, with the same
  shuffle and seed. Each region has its own exact null (only `permutations`,
  `shuffle` and `correction` apply). `_regions.tsv` gains `<metric>_p`,
  `_pfdr` and `_pfwer` columns, with matching statmaps. With an older rMVPA
  the analyses stay descriptive and `check` says what is missing. The studio
  marks, per parcel, the subjects whose own test survives FDR.
* Native-space projects: the studio finds the subject's fMRIPrep T1w and
  cortical surfaces (with bidser when installed, else by scanning
  `derivatives/*/sub-*/anat`) and draws native atlases on that subject's own
  surface, sampled at mid-thickness in the shared T1w space and tagged
  "sampled".

* The studio ships in the package (`inst/app`, one bundle built from
  `dsl/app`). `rmvpa_studio()` (`rmvpa studio`) serves it on 127.0.0.1 with a
  session token and Host/Origin checks; `rmvpa_export_studio()`
  (`rmvpa export <dir>`) writes it with the project's results to a folder
  that opens from disk with no R or server. Both read rmvpa_run() outputs only.
* The studio finds or fetches the images it draws, so none has to be named.
  Under the slices it offers the analysis mask, any T1w in the data space next
  to the data (fMRIPrep `anat/`), and for MNI spaces the TemplateFlow template,
  fetched when picked. Regional results go on fsaverage5: exactly for
  Schaefer-2018 atlases (CBIG's own annotation, recognised by file name or
  label table), approximately for other MNI-space atlases (sampled at
  mid-thickness, each hemisphere kept to its own parcels, and labelled
  approximate). Downloads happen only on the user's click (or
  `rmvpa_export_studio(fetch = TRUE)`) and are cached in
  `tools::R_user_dir("rmvpadsl", "cache")`.
* The studio has a whole-brain view, and httpuv, gifti and freesurferformats
  are now imported, so `rmvpa studio` works after a plain install.
* Searchlight permutation runs also write the pooled null as 1001 quantiles
  (`*_stat-<metric>_null.tsv`), so its distribution can be redrawn.

* First version of the rmvpa analysis language (DSL 0.1), moved out of rMVPA
  into its own package so that rMVPA's CRAN package stays free of it.
  Verbs: decode, transfer, rsa, rsa_vector, rsa_contrast, rsa_feature,
  rsa_era, rsa_connect, rsa_mediate, encode, pattern, manova and custom.
  Commands: new, add, explain, check, draft, run, status, lock and migrate.
