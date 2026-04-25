# Build inst/extdata/kriegeskorte92/rdms.rds
#
# Source: Kriegeskorte et al. (2008) Neuron 60:1126-1141, redistributed via
# the rsatoolbox demos/92imageData directory:
#   https://github.com/rsagroup/rsatoolbox/tree/main/demos/92imageData
#
# We bundle two derived files (brain RDMs, model RDMs) as a single compact RDS
# rather than the raw .mat supplemental package, which is ~3.5 MB.
#
# Run from the package root:
#   source("data-raw/kriegeskorte92.R")
#
# Requires: R.matlab, internet access (curl).

suppressPackageStartupMessages(library(R.matlab))

stopifnot(file.exists("DESCRIPTION"))

base_url <- "https://raw.githubusercontent.com/rsagroup/rsatoolbox/main/demos/92imageData"
tmp <- tempfile(fileext = "_kriegeskorte92")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

dl <- function(name) {
  dest <- file.path(tmp, name)
  utils::download.file(file.path(base_url, name), destfile = dest,
                       mode = "wb", quiet = TRUE)
  dest
}

bm <- readMat(dl("92_brainRDMs.mat"))
mm <- readMat(dl("92_modelRDMs.mat"))

brain_arr  <- bm$RDMs
brain_rdms <- list()
for (i in seq_len(dim(brain_arr)[3])) {
  for (j in seq_len(dim(brain_arr)[4])) {
    rdm   <- brain_arr[1, 1, i, j][[1]]
    label <- brain_arr[2, 1, i, j][[1]][1]
    nm    <- gsub("[^A-Za-z0-9]+", "_", label)
    brain_rdms[[nm]] <- as.matrix(rdm)
  }
}

model_arr  <- mm$Models
model_rdms <- list()
for (i in seq_len(dim(model_arr)[3])) {
  rdm  <- model_arr[1, 1, i][[1]]
  name <- model_arr[3, 1, i][[1]][1]
  model_rdms[[name]] <- as.matrix(rdm)
}

bundle <- list(
  brain    = brain_rdms,
  model    = model_rdms,
  source   = paste(
    "Kriegeskorte, Mur, Ruff, Kiani, Bodurka, Esteky, Tanaka, Bandettini (2008).",
    "Matching categorical object representations in inferior temporal cortex of man and monkey.",
    "Neuron 60:1126-1141. doi:10.1016/j.neuron.2008.10.043.",
    "Redistributed via rsatoolbox demos/92imageData (rsagroup/rsatoolbox)."
  ),
  n_items  = 92L,
  subjects = c("BE", "KO", "SN", "TI"),
  sessions = 1:2
)

dest_dir <- file.path("inst", "extdata", "kriegeskorte92")
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(bundle, file.path(dest_dir, "rdms.rds"), compress = "xz")

invisible(bundle)
