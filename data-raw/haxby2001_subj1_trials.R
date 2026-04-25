# Build inst/extdata/haxby2001_subj1/trials.rds
#
# Trial-level (per-TR) version of the Haxby 2001 Subject 1 ventral-temporal
# patterns. Sibling of patterns.rds (which stores per-(category, run) block
# means). The trial-level bundle is what you need for downstream multivariate
# noise normalisation (MVNN), since MVNN's noise covariance has to be
# estimated from within-block trial residuals -- block means alone do not
# carry that information.
#
# Source: same upstream archive as patterns.rds:
#   http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz
#
# Run from the package root:
#   source("data-raw/haxby2001_subj1_trials.R")
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
mat  <- do.call(cbind, lapply(vols, function(v) v[mask_idx]))   # V x T

# A "block" is a maximal run of consecutive same-category TRs in one chunk.
non_rest <- labels$labels != "rest"
block_id <- integer(nrow(labels))
cur  <- 0L
prev <- NA_character_
prev_chunk <- NA_integer_
for (i in seq_len(nrow(labels))) {
  if (labels$labels[i] == "rest") {
    prev <- NA_character_; prev_chunk <- NA_integer_; next
  }
  if (is.na(prev) ||
      labels$labels[i] != prev ||
      labels$chunks[i] != prev_chunk) {
    cur <- cur + 1L
  }
  block_id[i] <- cur
  prev <- labels$labels[i]
  prev_chunk <- labels$chunks[i]
}

keep <- non_rest
trial_X <- t(mat[, keep])
storage.mode(trial_X) <- "double"

bundle <- list(
  patterns   = trial_X,
  category   = factor(labels$labels[keep],
                      levels = setdiff(unique(labels$labels), "rest")),
  run        = as.integer(labels$chunks[keep]),
  block      = block_id[keep],
  mask_dim   = dim(mask),
  mask_idx   = mask_idx,
  mask_space = neuroim2::space(mask),
  source     = paste(
    "Haxby et al. (2001) Science 293:2425-2430.",
    "Subject 1 trial-level VT patterns from the PyMVPA-curated tutorial",
    "archive at http://data.pymvpa.org/datasets/haxby2001/subj1-2010.01.14.tar.gz.",
    "Each row is one task TR (rest TRs dropped); `block` is a within-run",
    "block index suitable for MVNN residual-covariance estimation."
  )
)

dest_dir <- file.path("inst", "extdata", "haxby2001_subj1")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, file.path(dest_dir, "trials.rds"), compress = "xz")

invisible(bundle)
