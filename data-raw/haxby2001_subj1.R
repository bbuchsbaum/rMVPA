# Build inst/extdata/haxby2001_subj1/patterns.rds
#
# Source: Haxby et al. (2001) Science 293:2425-2430. Subject 1 from the
# PyMVPA-curated tutorial archive:
#   http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz
#
# The upstream archive is ~300 MB (raw 4D BOLD + a VT mask + per-volume labels
# in attributes.txt). The bundle saved here is a small per-(category, run)
# mean-pattern matrix restricted to the published VT mask -- 8 categories x
# 12 runs x ~577 voxels, ~100 KB.
#
# Run from the package root:
#   source("data-raw/haxby2001_subj1.R")
#
# Requires: neuroim2, internet access (curl).

suppressPackageStartupMessages(library(neuroim2))

stopifnot(file.exists("DESCRIPTION"))

src_url <- "http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz"
tmp <- tempfile(fileext = "_haxby")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
archive <- file.path(tmp, "subj1.tar.gz")

utils::download.file(src_url, destfile = archive, mode = "wb", quiet = TRUE)
utils::untar(archive, exdir = tmp)

labels <- read.table(file.path(tmp, "subj1", "labels.txt"),
                     header = TRUE, stringsAsFactors = FALSE)
bold <- neuroim2::read_vec(file.path(tmp, "subj1", "bold.nii.gz"))
mask <- neuroim2::read_vol(file.path(tmp, "subj1", "mask4_vt.nii.gz"))

mask_idx <- which(mask > 0)
vols <- neuroim2::vols(bold)
mat  <- do.call(cbind, lapply(vols, function(v) v[mask_idx]))

non_rest <- setdiff(unique(labels$labels), "rest")
runs     <- sort(unique(labels$chunks))

n_obs <- length(non_rest) * length(runs)
X     <- matrix(NA_real_, nrow = n_obs, ncol = length(mask_idx))
cat_v <- character(n_obs)
rn_v  <- integer(n_obs)
k     <- 0L
for (c in non_rest) for (r in runs) {
  k <- k + 1L
  sel <- labels$labels == c & labels$chunks == r
  X[k, ]  <- rowMeans(mat[, sel, drop = FALSE])
  cat_v[k] <- c
  rn_v[k]  <- r
}
storage.mode(X) <- "double"

bundle <- list(
  patterns   = X,
  category   = factor(cat_v, levels = non_rest),
  run        = rn_v,
  mask_dim   = dim(mask),
  mask_idx   = mask_idx,
  mask_space = neuroim2::space(mask),
  source     = paste(
    "Haxby, Gobbini, Furey, Ishai, Schouten, Pietrini (2001).",
    "Distributed and overlapping representations of faces and objects in",
    "ventral temporal cortex. Science 293:2425-2430.",
    "Bundled here from the PyMVPA-curated subject-1 tutorial archive at",
    "http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz."
  )
)

dest_dir <- file.path("inst", "extdata", "haxby2001_subj1")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, file.path(dest_dir, "patterns.rds"), compress = "xz")

invisible(bundle)
