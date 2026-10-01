# rMVPA side of the benchmark harness. Run from the package root, after
# `Rscript tools/bench/make_data.R`:
#
#   Rscript tools/bench/run_rmvpa.R [reps] [groups]
#
# groups: comma-separated subset of sl,regional,rsa (default: all).
#
# Appends one JSON line per (scenario, method) to
# tools/bench/receipts/<date>-rmvpa.jsonl, in the same format as
# competitors/run_python.py. Runs are sequential. Engines are named
# explicitly, so a receipt states which code path it timed.

suppressMessages(pkgload::load_all(".", quiet = TRUE))
invisible(futile.logger::flog.threshold(futile.logger::INFO))

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args)) as.integer(args[[1]]) else 3L
groups <- if (length(args) >= 2) strsplit(args[[2]], ",")[[1]] else c("sl", "regional", "rsa")
data_dir <- file.path("tools", "bench", "data")
receipt_dir <- file.path("tools", "bench", "receipts")
dir.create(receipt_dir, recursive = TRUE, showWarnings = FALSE)
sha <- system("git rev-parse --short HEAD", intern = TRUE)
branch <- system("git rev-parse --abbrev-ref HEAD", intern = TRUE)

quiet <- function(code) {
  invisible(utils::capture.output(res <- suppressMessages(suppressWarnings(force(code)))))
  res
}

# `inner` repeats fn within each timed rep, for operations shorter than the
# timer resolution; times are reported per call.
#
# Every rep reseeds: rMVPA currently breaks predicted-class ties at random
# (release-plan finding C2), so without a fixed seed the output digest would
# depend on how much RNG earlier runs consumed.
time_reps <- function(fn, reps, inner = 1L) {
  set.seed(20261001)
  quiet(fn())  # warm-up
  out <- NULL
  times <- vapply(seq_len(reps), function(i) {
    set.seed(20261001)
    t0 <- proc.time()[["elapsed"]]
    for (k in seq_len(inner)) out <<- quiet(fn())
    (proc.time()[["elapsed"]] - t0) / inner
  }, numeric(1))
  list(times = times, out = out)
}

receipt <- function(scenario, method, times, n_units, unit, digest, engine = NA, notes = "") {
  list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    tool = "rMVPA",
    tool_versions = list(rMVPA = as.character(utils::packageVersion("rMVPA")),
                         neuroim2 = as.character(utils::packageVersion("neuroim2"))),
    r = R.version.string,
    blas = basename(extSoftVersion()[["BLAS"]]),
    threads = 1L,
    scenario = scenario,
    method = method,
    engine = engine,
    reps = length(times),
    median_s = stats::median(times),
    iqr_s = unname(diff(stats::quantile(times, c(0.25, 0.75)))),
    n_units = n_units,
    unit = unit,
    ms_per_unit = 1000 * stats::median(times) / n_units,
    digest = digest,
    notes = notes,
    rmvpa_git_sha = sha,
    rmvpa_branch = branch
  )
}

# ---- searchlight: synthetic 12^3, radius 3 mm (1 mm voxels) -----------------
bold <- neuroim2::read_vec(file.path(data_dir, "synth12_bold.nii.gz"))
mask <- neuroim2::read_vol(file.path(data_dir, "synth12_mask.nii.gz"))
sdes <- utils::read.csv(file.path(data_dir, "synth12_design.csv"))
sl_ds <- mvpa_dataset(bold, mask = neuroim2::LogicalNeuroVol(as.logical(mask), neuroim2::space(mask)))
sl_design <- mvpa_design(data.frame(label = factor(sdes$label), run = sdes$run),
                         y_train = ~ label, block_var = ~ run)
sl_cv <- blocked_cross_validation(sdes$run)
n_centres <- sum(as.logical(mask))

sl_methods <- list(
  corclass = list(model = "corclass", engine = "auto", note = "same model as nilearn corclass; best exact engine auto selects"),
  corclass_general = list(model = "corclass", engine = "legacy", note = "rMVPA only: general per-sphere path, for tracking"),
  gaussian_nb = list(model = "naive_bayes", engine = "auto", note = "same model as GaussianNB up to variance smoothing; best exact engine auto selects"),
  gaussian_nb_general = list(model = "naive_bayes", engine = "legacy", note = "rMVPA only: general per-sphere path, for tracking"),
  lda_shrinkage = list(model = "dual_lda", engine = "dual_lda_fast", note = "comparable, not identical: dual_lda gamma vs sklearn Ledoit-Wolf"),
  sda_notune = list(model = "sda_notune", engine = "auto", note = "rMVPA only; no sklearn equivalent; best exact engine auto selects"),
  sda_notune_general = list(model = "sda_notune", engine = "legacy", note = "rMVPA only: general per-sphere path, for tracking")
)

rows <- list()
for (nm in if ("sl" %in% groups) names(sl_methods) else character()) {
  m <- sl_methods[[nm]]
  mspec <- mvpa_model(load_model(m$model), sl_ds, sl_design, "classification", crossval = sl_cv)
  r <- time_reps(function() run_searchlight(mspec, radius = 3, method = "standard",
                                            engine = m$engine, backend = "default"), reps)
  acc <- as.numeric(neuroim2::values(r$out$results$Accuracy))[as.logical(mask)]
  rows[[length(rows) + 1]] <- receipt("sl_synth12_r3", nm, r$times, n_centres, "centre",
                                      list(mean_score = mean(acc, na.rm = TRUE)),
                                      engine = attr(r$out, "searchlight_engine"), notes = m$note)
}

# ---- regional: Haxby VT, leave-one-run-out ----------------------------------
bundle <- readRDS(system.file("extdata", "haxby2001_subj1", "patterns.rds", package = "rMVPA"))
mask_arr <- array(0L, bundle$mask_dim)
mask_arr[bundle$mask_idx] <- 1L
vec <- t(bundle$patterns)
storage.mode(vec) <- "double"
hx_ds <- mvpa_dataset(
  neuroim2::SparseNeuroVec(vec, neuroim2::add_dim(bundle$mask_space, ncol(vec)), mask = as.logical(mask_arr)),
  mask = neuroim2::LogicalNeuroVol(mask_arr, bundle$mask_space))
hx_design <- mvpa_design(data.frame(label = bundle$category, run = bundle$run),
                         y_train = ~ label, block_var = ~ run)
hx_cv <- blocked_cross_validation(bundle$run)
roi <- neuroim2::NeuroVol(mask_arr, bundle$mask_space)

for (nm in if ("regional" %in% groups) names(sl_methods) else character()) {
  m <- sl_methods[[nm]]
  mspec <- mvpa_model(load_model(m$model), hx_ds, hx_design, "classification", crossval = hx_cv)
  r <- time_reps(function() run_regional(mspec, roi), max(reps, 10L))
  rows[[length(rows) + 1]] <- receipt("regional_haxby_vt", nm, r$times, 1, "roi",
                                      list(accuracy = r$out$performance_table$Accuracy[[1]]),
                                      notes = m$note)
}

# ---- RSA: Haxby condition-mean correlation RDM and crossnobis ---------------
cats <- levels(bundle$category)
runs <- sort(unique(bundle$run))
if ("rsa" %in% groups) {
r <- time_reps(function() {
  means <- rowsum(bundle$patterns, bundle$category)[cats, ] / as.vector(table(bundle$category)[cats])
  d <- pairwise_dist(cordist(), means)
  d[lower.tri(d)]
}, max(reps, 20L), inner = 500L)
rows[[length(rows) + 1]] <- receipt("rsa_haxby", "rdm_correlation_condmeans", r$times, 1, "rdm",
                                    list(sum = sum(r$out)), notes = "8 condition means over 577 voxels")

r <- time_reps(function() {
  U <- array(NA_real_, dim = c(length(cats), ncol(bundle$patterns), length(runs)),
             dimnames = list(cats, NULL, NULL))
  for (m in seq_along(runs)) {
    sel <- bundle$run == runs[m]
    U[, , m] <- rowsum(bundle$patterns[sel, , drop = FALSE], bundle$category[sel])[cats, ] /
      as.vector(table(bundle$category[sel])[cats])
  }
  compute_crossnobis_distances_sl(U)
}, max(reps, 20L), inner = 100L)
rows[[length(rows) + 1]] <- receipt("rsa_haxby", "rdm_crossnobis_identity", r$times, 1, "rdm",
                                    list(sum = sum(r$out)), notes = "identity noise; runs as cv folds")
}

path <- file.path(receipt_dir, paste0(format(Sys.Date()), "-rmvpa.jsonl"))
for (row in rows) {
  cat(jsonlite::toJSON(row, auto_unbox = TRUE, digits = NA), "\n", file = path, append = TRUE, sep = "")
  cat(sprintf("%-18s %-28s rMVPA[%s] median=%.4fs  %.3f ms/%s  %s\n", row$scenario, row$method,
              row$engine %||% "-", row$median_s, row$ms_per_unit, row$unit,
              paste(names(row$digest), signif(unlist(row$digest), 6), collapse = " ")))
}
