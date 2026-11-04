# Synthetic BIDS-derivatives project for the rmvpa DSL tests.
#
# Two subjects, `nrun` runs of a task "objects" with 8 trials per run
# (2 categories x 4 stimuli, every stimulus once per run), a task "memory"
# with encoding/retrieval phases, a brain mask, a 2-region atlas and RDMs.

quiet_dsl <- function(expr) {
  futile.logger::flog.threshold(futile.logger::ERROR)
  on.exit(futile.logger::flog.threshold(futile.logger::INFO), add = TRUE)
  suppressWarnings(suppressMessages(utils::capture.output(res <- force(expr))))
  invisible(res)
}

dsl_skip_if_missing <- function() {
  testthat::skip_if_not_installed("yaml")
  testthat::skip_if_not_installed("jsonlite")
}

dsl_write_betas <- function(path_stem, design, space, signal_col, mask_arr, seed) {
  set.seed(seed)
  dims <- dim(mask_arr)
  n <- nrow(design)
  arr <- array(stats::rnorm(prod(dims) * n), c(dims, n))
  lev <- factor(design[[signal_col]])
  for (i in seq_len(n)) {
    s <- if (as.integer(lev[i]) == 1L) 0.8 else -0.8
    arr[2:4, 2:4, 2:4, i] <- arr[2:4, 2:4, 2:4, i] + s
  }
  vec <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(c(dims, n), spacing = c(2, 2, 2)))
  neuroim2::write_vec(vec, paste0(path_stem, ".nii.gz"))
  utils::write.table(design, paste0(path_stem, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

dsl_fixture <- function(dir = tempfile("rmvpa-dsl-"), subjects = c("01", "02"), nrun = 4L) {
  dims <- c(6L, 6L, 6L)
  space <- neuroim2::NeuroSpace(dims, spacing = c(2, 2, 2))
  mask_arr <- array(1, dims)
  bids <- file.path(dir, "derivatives", "lss")
  stims <- c("cat1", "cat2", "dog1", "dog2", "car1", "car2", "bus1", "bus2")
  category <- c("animal", "animal", "animal", "animal", "vehicle", "vehicle", "vehicle", "vehicle")
  for (si in seq_along(subjects)) {
    s <- subjects[si]
    fdir <- file.path(bids, paste0("sub-", s), "func")
    dir.create(fdir, recursive = TRUE, showWarnings = FALSE)
    neuroim2::write_vol(neuroim2::NeuroVol(mask_arr, space),
                        file.path(fdir, sprintf("sub-%s_space-MNI_desc-brain_mask.nii.gz", s)))
    for (r in seq_len(nrun)) {
      d <- data.frame(stimulus = stims, category = category,
                      response = rep(c("hit", "hit", "hit", "miss"), 2), stringsAsFactors = FALSE)
      dsl_write_betas(file.path(fdir, sprintf("sub-%s_task-objects_run-%d_space-MNI_desc-lss_betas", s, r)),
                      d, space, "category", mask_arr, seed = 100 * si + r)
      for (ph in c("encoding", "retrieval")) {
        dm <- data.frame(stimulus = stims, category = category, phase = ph, stringsAsFactors = FALSE)
        dsl_write_betas(file.path(fdir, sprintf("sub-%s_task-memory_acq-%s_run-%d_space-MNI_desc-lss_betas",
                                                s, ph, r)),
                        dm, space, "category", mask_arr, seed = 1000 * si + 10 * r + nchar(ph))
      }
    }
  }
  # Atlas: two regions.
  atlas <- array(0L, dims); atlas[1:3, , ] <- 1L; atlas[4:6, , ] <- 2L
  dir.create(file.path(dir, "atlas"), showWarnings = FALSE)
  neuroim2::write_vol(neuroim2::NeuroVol(atlas, space), file.path(dir, "atlas", "space-MNI_atlas-two_dseg.nii.gz"))
  utils::write.table(data.frame(index = 1:2, name = c("left", "right")),
                     file.path(dir, "atlas", "space-MNI_atlas-two_dseg.tsv"), sep = "\t", quote = FALSE,
                     row.names = FALSE)
  # RDMs over stimuli.
  dir.create(file.path(dir, "rdms"), showWarnings = FALSE)
  cat_rdm <- outer(category, category, function(a, b) as.numeric(a != b))
  dimnames(cat_rdm) <- list(stims, stims)
  set.seed(7)
  noise <- as.matrix(stats::dist(stats::rnorm(8))); dimnames(noise) <- list(stims, stims)
  write_rdm <- function(m, f) utils::write.csv(data.frame(item = rownames(m), m, check.names = FALSE),
                                               f, row.names = FALSE)
  write_rdm(cat_rdm, file.path(dir, "rdms", "category.csv"))
  write_rdm(noise, file.path(dir, "rdms", "noise.csv"))
  dir
}

dsl_write_project <- function(dir, extra_sets = TRUE, assert = "classes_in_every_fold") {
  lines <- c(
    'rmvpa: "0.1"',
    "data:",
    "  bids: derivatives/lss",
    "  space: MNI",
    "  mask: brain",
    "  sets:",
    "    objects: { task: objects, desc: lss }"
  )
  if (extra_sets) {
    lines <- c(lines,
               "    encoding: { task: memory, acq: encoding }",
               "    retrieval: { task: memory, acq: retrieval }")
  }
  lines <- c(lines, "folds:", "  by: run", if (length(assert)) sprintf("  assert: [%s]", paste(assert, collapse = ", ")),
             "run:", "  seed: 42")
  writeLines(lines, file.path(dir, "rmvpa.yml"))
  dir.create(file.path(dir, "analyses"), showWarnings = FALSE)
  invisible(dir)
}

dsl_write_analysis <- function(dir, name, lines) {
  writeLines(c("extends: ../rmvpa.yml", lines), file.path(dir, "analyses", paste0(name, ".yml")))
  invisible(file.path(dir, "analyses", paste0(name, ".yml")))
}
