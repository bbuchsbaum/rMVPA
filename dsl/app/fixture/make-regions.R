# Regional example for the studio: transfer_memory on the Schaefer-200
# (7-network) atlas, shown on slices and on the fsaverage5 cortical surface.
#
# The atlas and meshes are real (CBIG Schaefer 2018; TemplateFlow fsaverage5
# pial + sulc). Accuracies are SYNTHETIC example data.
#
#   Rscript dsl/app/fixture/make-regions.R <dir with downloaded files>
#
# Expected files in <dir>:
#   Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz
#   {lh,rh}.Schaefer2018_200Parcels_7Networks_order.annot      (fsaverage5)
#   tpl-fsaverage_hemi-{L,R}_den-10k_{pial.surf,sulc.shape}.gii
suppressMessages({ library(neuroim2) })
set.seed(4511)
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else "."
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
out <- normalizePath(file.path(here, "..", "public", "example"), mustWork = FALSE)
b64 <- function(x, size = 4, int = FALSE) jsonlite::base64_enc(writeBin(if (int) as.integer(x) else as.numeric(x), raw(), size = size))

# ---- volume: Schaefer labels resampled (nearest) onto the studio template ----
tpl <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
sch <- read_vol(file.path(src, "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"))
d <- dim(tpl)
g <- as.matrix(expand.grid(1:d[1], 1:d[2], 1:d[3]))
# neuroim2's mni_downsampled header records 1 mm voxels, but the 48 x 57 x 48
# grid is the MNI brain at 4 mm; use the 4 mm affine with the header origin.
origin <- trans(space(tpl))[1:3, 4]
world <- sweep((g - 1) * 4, 2, origin, "+")
sg <- round(coord_to_grid(space(sch), world))
ok <- sg[, 1] >= 1 & sg[, 2] >= 1 & sg[, 3] >= 1 & sg[, 1] <= dim(sch)[1] & sg[, 2] <= dim(sch)[2] & sg[, 3] <= dim(sch)[3]
lab <- integer(nrow(g)); lab[ok] <- as.integer(as.array(sch)[sg[ok, , drop = FALSE]])
cat("voxels labelled:", sum(lab > 0), "parcels present:", length(unique(lab[lab > 0])), "\n")

# ---- regions from the annotation colour tables ----
annot <- lapply(c(lh = "lh", rh = "rh"), function(h)
  freesurferformats::read.fs.annot(file.path(src, sprintf("%s.Schaefer2018_200Parcels_7Networks_order.annot", h))))
names_lh <- setdiff(annot$lh$colortable$struct_names, c("Background+FreeSurfer_Defined_Medial_Wall", "Unknown", "???"))
names_rh <- setdiff(annot$rh$colortable$struct_names, c("Background+FreeSurfer_Defined_Medial_Wall", "Unknown", "???"))
full <- c(names_lh, names_rh)
stopifnot(length(full) == 200)
net <- sub("^7Networks_[LR]H_([A-Za-z]+)_.*$", "\\1", full)
regions <- data.frame(id = 1:200, name = full, hemisphere = ifelse(grepl("_LH_", full), "left", "right"),
                      network = net, label = sub("^7Networks_[LR]H_", "", full), stringsAsFactors = FALSE)

# ---- synthetic transfer accuracies (16 subjects, sub-16 included: all cells done) ----
chance <- 0.25
eff <- c(Vis = 0.11, DorsAttn = 0.045, SomMot = 0.01, SalVentAttn = 0.012, Limbic = 0, Cont = 0.015, Default = 0.008)
centroid <- t(sapply(1:200, function(r) { w <- which(lab == r); if (!length(w)) c(NA, NA, NA) else colMeans(world[w, , drop = FALSE]) }))
post <- ifelse(is.na(centroid[, 2]), 0, pmax(0, (-centroid[, 2] - 40) / 60))   # posterior gradient within networks
mu <- chance + eff[regions$network] * (0.6 + 0.8 * post) + rnorm(200, 0, 0.006)
gain <- rnorm(16, 1, 0.3)
vals <- sapply(1:16, function(s) round(chance + (mu - chance) * gain[s] + rnorm(200, 0, 0.03), 4))
region_rows <- lapply(1:200, function(r) list(index = r, name = regions$label[r], full = regions$name[r],
                                              hemisphere = regions$hemisphere[r], network = regions$network[r],
                                              values = vals[r, ]))

# ---- surface: fsaverage5 pial, sulc, and per-vertex parcel ids ----
hemi <- function(H, h) {
  gs <- gifti::read_gifti(file.path(src, sprintf("tpl-fsaverage_hemi-%s_den-10k_pial.surf.gii", H)))
  pos <- gs$data$pointset; tri <- gs$data$triangle
  sulc <- gifti::read_gifti(file.path(src, sprintf("tpl-fsaverage_hemi-%s_den-10k_sulc.shape.gii", H)))$data[[1]]
  a <- annot[[h]]
  ids <- match(a$label_names, full); ids[is.na(ids)] <- 0L
  stopifnot(nrow(pos) == length(ids))
  list(n = nrow(pos), positions = b64(t(pos)), faces = b64(t(tri), int = TRUE),
       sulc = b64(as.numeric(sulc)), labels = b64(ids, int = TRUE))
}
# atlas_md5 ties the vertex labels to the volumetric atlas file whose parcel
# ids they use: rmvpa_studio(surfaces = ) draws regional results on this mesh
# only for analyses run on that atlas.
surface <- list(template = "fsaverage5", mesh = "pial", source = "TemplateFlow tpl-fsaverage den-10k; CBIG Schaefer2018 fsaverage5 annot",
                atlas_md5 = unname(tools::md5sum(file.path(src, "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"))),
                left = hemi("L", "lh"), right = hemi("R", "rh"))

dir.create(out, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(list(
  analysis = "transfer_memory", metric = "Accuracy", chance = chance,
  chance_source = "1 / 4 classes; the retrieval set is balanced (36 trials per class per subject)",
  subjects = sprintf("%02d", 1:16),
  atlas = list(name = "Schaefer2018 200 parcels, 7 networks", file = "atlas/space-MNI152NLin6Asym_atlas-Schaefer2018_desc-200Parcels7Networks_dseg.nii.gz",
               md5 = unname(tools::md5sum(file.path(src, "Schaefer2018_200Parcels_7Networks_order_FSLMNI152_2mm.nii.gz"))),
               surface = "fsaverage5"),
  labels = b64(lab, int = TRUE), regions = region_rows,
  null = list(kind = "none", description = "none: descriptive only (no permutation path in rMVPA for this analysis)")
), file.path(out, "results-transfer_memory.json"), auto_unbox = TRUE, digits = 5)
jsonlite::write_json(surface, file.path(out, "surface-fsaverage5.json"), auto_unbox = TRUE)
cat("Wrote", out, "\n")
