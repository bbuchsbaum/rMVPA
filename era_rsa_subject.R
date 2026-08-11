#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(neuroim2); library(rMVPA); library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
sid <- if (length(args)) args[[1]] else "1002"
root <- "/scratch/brad/dsets/pcmri"
func <- file.path(root, "derivatives/lss_trialwise", paste0("sub-", sid), "func")
out <- file.path(root, "derivatives/era_rsa", paste0("sub-", sid))
cores <- as.integer(Sys.getenv("ERA_RSA_CORES", "1"))
radius <- as.numeric(Sys.getenv("ERA_RSA_RADIUS_MM", "2"))
batch_size <- as.integer(Sys.getenv("ERA_RSA_BATCH_SIZE", "64"))
if (is.na(batch_size) || batch_size < 1L) stop("ERA_RSA_BATCH_SIZE must be a positive integer")
options(future.globals.maxSize = 2 * 1024^3)
mask_file <- file.path(root, "derivatives/glm_recog", paste0("sub-", sid), paste0("sub-", sid, "_task-recog_desc-brain_mask.nii.gz"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stem <- function(task, desc) file.path(func, sprintf("sub-%s_task-%s_desc-%s_bold", sid, task, desc))
enc_stem <- stem("encode", "lssEncoding")
ret_stem <- stem("recog", "lssRetrieval")
enc_file <- paste0(enc_stem, ".nii.gz")
ret_file <- paste0(ret_stem, ".nii.gz")
enc_meta <- read_tsv(paste0(enc_stem, "_trials.tsv"), show_col_types = FALSE) |> mutate(ImageNumber = as.character(ImageNumber), enc_row = row_number())
ret_meta <- read_tsv(paste0(ret_stem, "_trials.tsv"), show_col_types = FALSE) |> mutate(ImageNumber = as.character(ImageNumber), Saliency = as.numeric(Saliency), ret_row = row_number())
enc_items <- enc_meta |> group_by(ImageNumber) |> summarise(enc_rows = list(enc_row), n_enc = n(), .groups = "drop")
match_tbl <- ret_meta |> filter(Condition == "old") |> left_join(enc_items, by = "ImageNumber") |> filter(!is.na(n_enc))
if (!nrow(match_tbl)) stop("No old-item matches for ", sid)
manifest_tbl <- match_tbl |> mutate(enc_rows = vapply(enc_rows, paste, character(1), collapse = ","))
write_tsv(manifest_tbl |> select(subject, trial_index, ImageNumber, Saliency, Run, enc_rows, n_enc, eye_sim_diff, sim_group_diff), file.path(out, "matching_manifest.tsv"))

stopifnot(file.exists(mask_file), file.exists(enc_file), file.exists(ret_file))
mask_vol <- read_vol(mask_file)
mask_values <- as.numeric(values(mask_vol))
mask_arr <- is.finite(mask_values) & mask_values > 0
dim(mask_arr) <- dim(mask_vol)
mask <- DenseNeuroVol(mask_arr, space(mask_vol), label = "era_rsa_brain_mask")

ret_rows <- match(match_tbl$trial_index, ret_meta$trial_index)
if (anyNA(ret_rows)) stop("Retrieval trial indices failed to map")
enc_rows_needed <- sort(unique(unlist(match_tbl$enc_rows, use.names = FALSE)))
match_tbl$enc_rows <- lapply(match_tbl$enc_rows, function(rows) {
  loaded_rows <- match(rows, enc_rows_needed)
  if (anyNA(loaded_rows)) stop("Encoding trial indices failed to map")
  loaded_rows
})

# Load only analysis-relevant volumes and only masked voxels. Encoding and
# retrieval remain separate throughout; no dense 4D concatenation is needed.
enc <- read_vec(enc_file, indices = enc_rows_needed, mask = mask_arr, mode = "bigvec")
ret <- read_vec(ret_file, indices = ret_rows, mask = mask_arr, mode = "bigvec")
dset <- mvpa_dataset(enc, test_data = ret, mask = mask)

n_encoding <- nrow(enc_meta)
n_encoding_selected <- length(enc_rows_needed)
n_retrieval <- nrow(ret_meta)
n_retrieval_selected <- length(ret_rows)
n_old_matches <- nrow(match_tbl)

# This remains an ordinary arbitrary custom callback. rMVPA supplies training
# data as `sl_data`, the corresponding retrieval sphere as
# `sl_info$test_data`, and the caller's item metadata as `sl_info$user_data`.
make_effects <- function(sl_data, sl_info) {
  encoding <- as.matrix(sl_data)
  retrieval <- sl_info$test_data
  mt <- sl_info$user_data
  if (is.null(retrieval)) stop("ERA-RSA requires a separate retrieval test set")
  if (nrow(retrieval) != nrow(mt)) stop("Retrieval rows and item metadata differ")

  enough_voxels <- ncol(encoding) >= 3L
  scores <- vapply(seq_len(nrow(mt)), function(i) {
    if (!enough_voxels) return(NA_real_)
    e <- colMeans(encoding[mt$enc_rows[[i]], , drop = FALSE], na.rm = TRUE)
    r <- retrieval[i, ]
    finite <- is.finite(e) & is.finite(r)
    if (sum(finite) < 3L || sd(e[finite]) == 0 || sd(r[finite]) == 0) {
      NA_real_
    } else {
      cor(e[finite], r[finite], method = "spearman")
    }
  }, numeric(1))

  valid_scores <- is.finite(scores)
  saliency_ok <- valid_scores & is.finite(mt$Saliency)
  safe_cor <- function(x, y, ok, method = "spearman", min_n = 4L) {
    if (sum(ok) < min_n || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
    cor(x[ok], y[ok], method = method)
  }
  ans <- list(
    n_vox = ncol(encoding),
    n_all = sum(valid_scores),
    mean_reactivation = if (any(valid_scores)) mean(scores[valid_scores]) else NA_real_,
    saliency_cor = safe_cor(scores, mt$Saliency, saliency_ok),
    saliency_n = sum(saliency_ok)
  )

  for (nm in c("eye_sim_diff", "sim_group_diff")) {
    z <- valid_scores & is.finite(mt[[nm]])
    ans[[paste0(nm, "_beta")]] <- if (sum(z) >= 4L) {
      unname(coef(lm(scores[z] ~ mt[[nm]][z]))[2])
    } else {
      NA_real_
    }
    ans[[paste0(nm, "_n")]] <- sum(z)
  }
  ans
}

sl <- run_custom_searchlight(
  dset,
  make_effects,
  radius = radius,
  method = "standard",
  user_data = match_tbl,
  batch_size = batch_size,
  .cores = cores,
  .verbose = TRUE
)
save_results(sl, out, level = "standard", stack = "auto", fname = "era_rsa_searchlight.nii.gz", overwrite = TRUE, quiet = FALSE)
# Keep the searchlight object as the reproducibility artifact; save_results()
# writes the metric NIfTIs above.
saveRDS(sl, file.path(out, "searchlight_result.rds"))
writeLines(c(paste0("subject=", sid), paste0("rMVPA=", packageVersion("rMVPA")), paste0("cores=", cores), paste0("batch_size=", batch_size), paste0("radius_mm=", radius), paste0("n_encoding=", n_encoding), paste0("n_encoding_selected=", n_encoding_selected), paste0("n_retrieval=", n_retrieval), paste0("n_retrieval_selected=", n_retrieval_selected), paste0("n_old_matches=", n_old_matches), "subject_level_inference=effect_estimation_only", "second_level_inference=deferred"), file.path(out, "era_rsa_summary.txt"))
cat("ERA-RSA production complete for sub-", sid, "\n", sep = "")
