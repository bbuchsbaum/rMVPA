# Generate benchmark inputs shared by the rMVPA and competitor runners.
# Run from the package root:  Rscript tools/bench/make_data.R
# Writes tools/bench/data/ (gitignored): a synthetic 4D NIfTI volume, its mask,
# trial labels/runs, and the Haxby VT block patterns as CSV.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
out <- file.path("tools", "bench", "data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

set.seed(20261001)
D <- c(12, 12, 12)
ds <- gen_sample_dataset(D = D, nobs = 100, nlevels = 4, blocks = 5)
neuroim2::write_vec(ds$dataset$train_data, file.path(out, "synth12_bold.nii.gz"))
mask <- neuroim2::NeuroVol(as.numeric(as.logical(ds$dataset$mask)), neuroim2::space(ds$dataset$mask))
neuroim2::write_vol(mask, file.path(out, "synth12_mask.nii.gz"))
utils::write.csv(data.frame(label = as.character(ds$design$y_train),
                            run = ds$design$block_var),
                 file.path(out, "synth12_design.csv"), row.names = FALSE)

# RSA searchlight scenario: the first 20 synthetic observations act as 20
# conditions; a fixed random model RDM over them (190 pairs, lower triangle
# column-major = upper triangle row-major).
set.seed(20261002)
utils::write.csv(data.frame(model = stats::runif(190)), file.path(out, "rsa_model_rdm.csv"), row.names = FALSE)

bundle <- readRDS(system.file("extdata", "haxby2001_subj1", "patterns.rds", package = "rMVPA"))
utils::write.csv(bundle$patterns, file.path(out, "haxby_patterns.csv"), row.names = FALSE)
utils::write.csv(data.frame(label = as.character(bundle$category), run = bundle$run),
                 file.path(out, "haxby_design.csv"), row.names = FALSE)
cat("wrote", out, "\n")
