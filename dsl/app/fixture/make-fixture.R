# Builds the studio's example project: dsl/app/public/example/*.json.
#
# The project files, `explain` output and `check` findings are real: they come
# from rmvpadsl run on a demo project whose design tables and mask exist on
# disk. The result maps and per-subject numbers are SYNTHETIC (example data for
# designing the interface), drawn on neuroim2's downsampled MNI template.
#
#   Rscript dsl/app/fixture/make-fixture.R
suppressMessages({ library(rmvpadsl); library(neuroim2) })
set.seed(20261104)
`%||%` <- function(a, b) if (is.null(a)) b else a
out <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))), "..", "public", "example"), mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
b64 <- function(x, what = "double", size = 4) jsonlite::base64_enc(writeBin(x, raw(), size = size))
save_json <- function(x, f) jsonlite::write_json(x, file.path(out, f), auto_unbox = TRUE, digits = 5, null = "null")

# ---- template and mask -------------------------------------------------------
tpl <- read_vol(system.file("extdata", "mni_downsampled.nii.gz", package = "neuroim2"))
dims <- dim(tpl); arr <- as.array(tpl)
mask <- arr > quantile(arr[arr > 0], 0.25)
tpl8 <- as.integer(pmin(255, round(255 * arr / quantile(arr, 0.995))))

# ---- demo project on disk (designs + mask real; betas are placeholders) -------
root <- file.path(tempdir(), "objects-study"); unlink(root, recursive = TRUE)
subs <- sprintf("%02d", 1:16)
cats <- c("face", "animal", "vehicle", "tool")
stims <- paste0(rep(cats, each = 6), rep(1:6, 4))
sp <- NeuroSpace(dims, spacing = c(4, 4, 4))
for (s in subs) {
  fd <- file.path(root, "derivatives", "lss", paste0("sub-", s), "func"); dir.create(fd, recursive = TRUE)
  write_vol(NeuroVol(array(as.numeric(mask), dims), sp), file.path(fd, sprintf("sub-%s_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz", s)))
  for (r in 1:6) {
    d <- data.frame(stimulus = sample(stims), stringsAsFactors = FALSE)
    d$category <- sub("[0-9]+$", "", d$stimulus)
    d$response <- ifelse(runif(nrow(d)) < 0.9, "hit", "miss")
    d$rt <- round(rlnorm(nrow(d), log(0.65), 0.25), 3)
    stem <- sprintf("sub-%s_task-objects_run-%d_space-MNI152NLin2009cAsym_desc-lss_betas", s, r)
    write.table(d, file.path(fd, paste0(stem, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
    file.create(file.path(fd, paste0(stem, ".nii.gz")))
    for (ph in c("encoding", "retrieval")) {
      dm <- data.frame(stimulus = sample(stims), stringsAsFactors = FALSE); dm$category <- sub("[0-9]+$", "", dm$stimulus)
      stem2 <- sprintf("sub-%s_task-memory_acq-%s_run-%d_space-MNI152NLin2009cAsym_desc-lss_betas", s, ph, r)
      write.table(dm, file.path(fd, paste0(stem2, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
      file.create(file.path(fd, paste0(stem2, ".nii.gz")))
    }
  }
}
dir.create(file.path(root, "atlas")); dir.create(file.path(root, "rdms")); dir.create(file.path(root, "analyses"))
# Atlas: 24 parcels by octant-ish slabs of the mask.
lab <- array(0L, dims)
idx <- which(mask, arr.ind = TRUE)
cx <- cut(idx[, 1], c(-Inf, dims[1] / 2, Inf), labels = FALSE); cy <- cut(idx[, 2], quantile(idx[, 2], 0:6 / 6), include.lowest = TRUE, labels = FALSE)
cz <- cut(idx[, 3], quantile(idx[, 3], 0:2 / 2), include.lowest = TRUE, labels = FALSE)
lab[idx] <- (cx - 1L) * 12L + (cy - 1L) * 2L + cz
write_vol(NeuroVol(array(as.numeric(lab), dims), sp), file.path(root, "atlas", "space-MNI152NLin2009cAsym_atlas-demo_dseg.nii.gz"))
lobes <- c("Occipital", "Occipital", "Temporal", "Parietal", "Frontal", "Frontal")
part <- c("posterior", "anterior", "", "", "posterior", "anterior")
regions <- data.frame(index = 1:24, name = trimws(gsub("  +", " ", sprintf("%s %s %s %s", rep(c("LH", "RH"), each = 12),
  rep(rep(lobes, each = 2), 2), rep(rep(part, each = 2), 2), rep(c("inferior", "superior"), 12)))))
write.table(regions, file.path(root, "atlas", "space-MNI152NLin2009cAsym_atlas-demo_dseg.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
anim <- outer(stims, stims, function(a, b) as.numeric(grepl("^(face|animal)", a) != grepl("^(face|animal)", b)))
catr <- outer(stims, stims, function(a, b) as.numeric(sub("[0-9]+$", "", a) != sub("[0-9]+$", "", b)))
wr <- function(m, f) write.csv(data.frame(item = stims, m, check.names = FALSE), file.path(root, "rdms", f), row.names = FALSE)
dimnames(anim) <- dimnames(catr) <- list(stims, stims); wr(anim, "animacy.csv"); wr(catr, "category.csv")
tm <- as.matrix(dist(seq_along(stims))) / length(stims); dimnames(tm) <- list(stims, stims); wr(tm, "time.csv")

writeLines(c('rmvpa: "0.1"', "data:", "  bids: derivatives/lss", "  space: MNI152NLin2009cAsym", "  mask: brain", "  sets:",
  "    objects:", "      task: objects", "      keep: { response: hit }",
  "    encoding:  { task: memory, acq: encoding }", "    retrieval: { task: memory, acq: retrieval }",
  "folds:", "  by: run", "  assert: [classes_in_every_fold]", "run:", "  workers: 8", "  seed: 20261104"),
  file.path(root, "rmvpa.yml"))
A <- list(
  decode_category = c("extends: ../rmvpa.yml", "decode:", "  set: objects", "  target: category", "scope:", "  searchlight: { radius: 6 }"),
  rsa_animacy = c("extends: ../rmvpa.yml", "rsa:", "  set: objects", "  item: stimulus", "  predictors:", "    animacy: rdms/animacy.csv",
                  "    category: rdms/category.csv", "  nuisance:", "    time: rdms/time.csv", "  method: lm", "scope:", "  searchlight: { radius: 8 }",
                  "inference:", "  permutations: 1000", "folds:", "  assert: ~"),
  transfer_memory = c("extends: ../rmvpa.yml", "transfer:", "  train_set: encoding", "  test_set: retrieval", "  target: category",
                      "scope:", "  regions:", "    atlas: atlas/space-MNI152NLin2009cAsym_atlas-demo_dseg.nii.gz",
                      "    labels: atlas/space-MNI152NLin2009cAsym_atlas-demo_dseg.tsv"),
  era_items = c("extends: ../rmvpa.yml", "rsa_era:", "  train_set: encoding", "  test_set: retrieval", "  link: stimulus",
                "folds:", "  assert: [no_item_across_folds]",
                "scope:", "  regions:", "    atlas: atlas/space-MNI152NLin2009cAsym_atlas-demo_dseg.nii.gz")
)
for (n in names(A)) writeLines(A[[n]], file.path(root, "analyses", paste0(n, ".yml")))

proj <- rmvpa_project(root)
chk <- as.data.frame(unclass(rmvpa_check(proj, quiet = TRUE)))

explain_rows <- function(spec) {
  rows <- list()
  walk <- function(x, path, depth) {
    for (k in names(x)) {
      p <- if (nzchar(path)) paste0(path, ".", k) else k
      v <- x[[k]]
      if (rmvpadsl:::.dsl_is_map(v) && length(v) && !identical(names(v), "in")) {
        rows[[length(rows) + 1L]] <<- list(path = p, key = k, depth = depth, value = NULL, origin = NULL)
        walk(v, p, depth + 1L)
      } else {
        o <- spec$origins[p]; o <- if (is.na(o)) NULL else unname(o)
        rows[[length(rows) + 1L]] <<- list(path = p, key = k, depth = depth, value = rmvpadsl:::.dsl_inline(v), origin = o)
      }
    }
  }
  cfg <- spec$config[intersect(c(spec$verb, "folds", "scope", "inference", "run"), names(spec$config))]
  walk(cfg, "", 0L)
  rows
}

# ---- synthetic results -------------------------------------------------------
blob <- function(c, s) { g <- expand.grid(1:dims[1], 1:dims[2], 1:dims[3]); array(exp(-((g[,1]-c[1])^2 + (g[,2]-c[2])^2 + (g[,3]-c[3])^2) / (2 * s^2)), dims) }
mid <- dims[1] / 2
vot <- blob(c(mid - 9, 16, 14), 3.2) + blob(c(mid + 9, 16, 14), 3.2) + 0.6 * blob(c(mid - 7, 9, 20), 3) + 0.6 * blob(c(mid + 7, 9, 20), 3) +
  0.35 * blob(c(mid, 8, 24), 4) + 0.25 * blob(c(mid - 14, 36, 30), 3)
smooth_noise <- function(sd) { n <- array(rnorm(prod(dims), 0, sd), dims); for (k in 1:2) n <- (n + aperm(apply(n, c(2, 3), function(v) stats::filter(v, rep(1/3, 3), circular = TRUE)), c(1, 2, 3))) / 2; n }
chance <- 0.25
grp <- chance + 0.14 * vot + smooth_noise(0.012); grp[!mask] <- NA
peak_idx <- which.max(grp)

nperm <- 100L
n_centres <- sum(mask)
n_null <- nperm * round(0.1 * n_centres)  # pooled null: permutations x subsampled centres
null_all <- rnorm(n_null, chance, 0.028) + rgamma(n_null, 2, 400)
null_sorted <- sort(null_all)
obs <- grp[mask]
pv <- (1 + (n_null - findInterval(obs, null_sorted))) / (n_null + 1)
edges <- seq(0.10, 0.50, by = 0.004)
null_hist <- as.integer(table(cut(pmin(pmax(null_all, 0.1), 0.4999), edges, include.lowest = TRUE)))
null_vals <- round(sample(null_all, 3000), 4)
q <- p.adjust(pv, "BH")
fdr_thr <- if (any(q < 0.05)) min(obs[q < 0.05]) else NA
pmap <- array(NA_real_, dims); pmap[mask] <- pv
# Per-subject maps at the peak slice; leave-one-subject-out ROI values.
pk <- arrayInd(which.max(grp), dims); kz <- pk[3]
gains <- rnorm(length(subs), 1, 0.35)
prev <- array(0L, dims); prev[!mask] <- NA   # subjects significant, done cells only (stale sub-07 excluded)
prev_all <- prev                                # including the stale subject
per <- vector("list", length(subs))
for (i in seq_len(15)) {
  m <- chance + 0.14 * gains[i] * vot + array(rnorm(prod(dims), 0, 0.035), dims)
  ov <- m[mask]
  nl <- sort(rnorm(n_null, chance, 0.03) + rgamma(n_null, 2, 400))
  p_i <- (1 + (n_null - findInterval(ov, nl))) / (n_null + 1)
  q_i <- p.adjust(p_i, "BH")
  sig <- array(FALSE, dims); sig[mask] <- q_i < 0.05
  prev_all[sig] <- prev_all[sig] + 1L
  if (i != 7) prev[sig] <- prev[sig] + 1L
  sl <- m[, , kz]; sl[!mask[, , kz]] <- NaN
  per[[i]] <- list(
    fdr = list(threshold = if (any(q_i < 0.05)) round(min(ov[q_i < 0.05]), 4) else NULL, surviving = sum(q_i < 0.05)),
    centres_hist = as.integer(table(cut(pmin(pmax(ov, 0.1), 0.4999), edges, include.lowest = TRUE))),
    null_hist = as.integer(table(cut(pmin(pmax(nl, 0.1), 0.4999), edges, include.lowest = TRUE))),
    slice = jsonlite::base64_enc(writeBin(as.numeric(sl), raw(), size = 4)))
}
subj_slice <- lapply(per, function(x) x$slice)
subj <- lapply(seq_along(subs), function(i) {
  list(id = subs[i], loso_roi = round(chance + 0.075 * gains[i] + rnorm(1, 0, 0.018), 4),
       status = "done", seconds = round(runif(1, 380, 560)), n_trials = 144L, kept = 144L - rbinom(1, 144, 0.1),
       slice = subj_slice[[i]], fdr = per[[i]]$fdr, centres_hist = per[[i]]$centres_hist, null_hist = per[[i]]$null_hist)
})
subj[[7]]$status <- "stale"; subj[[16]]$status <- "missing"; subj[[16]]$loso_roi <- NULL

rsa_t <- function(w, s) { m <- w * vot + smooth_noise(s); m[!mask] <- NA; m }
anim_t <- rsa_t(5.2, 0.5); cat_t <- rsa_t(3.0, 0.5); time_t <- rsa_t(0.4, 0.6)

reg_vals <- function(eff) {
  sapply(seq_len(24), function(r) { v <- vot[lab == r]; if (!length(v)) 0 else mean(v) }) * eff
}
xfer <- reg_vals(0.9)
xfer_tab <- lapply(1:24, function(r) list(index = r, name = regions$name[r],
  values = round(chance + xfer[r] + rnorm(16, 0, 0.035), 4)))
era <- reg_vals(0.5)
era_tab <- lapply(1:24, function(r) list(index = r, name = regions$name[r], values = round(era[r] + rnorm(15, 0, 0.03), 4)))

status_of <- function(name) {
  st <- rep("done", length(subs))
  if (name == "decode_category") st[c(7, 16)] <- c("stale", "missing")
  if (name == "rsa_animacy") st[11:16] <- "missing"
  if (name == "era_items") st[] <- "blocked"
  st
}

analyses <- lapply(names(A), function(n) {
  spec <- rmvpa_analysis(proj, n)
  list(name = n, verb = spec$verb, family = spec$family, scope = spec$scope_kind,
       scope_detail = spec$config$scope[[1]], file = paste(A[[n]], collapse = "\n"),
       explain = explain_rows(spec), notes = unlist(spec$notes),
       inference = list(kind = spec$inference_support$kind, description = spec$inference_support$description,
                        test = spec$config$inference$test, permutations = spec$config$inference$permutations %||% NULL),
       lock = c(decode_category = "match", rsa_animacy = "mismatch", transfer_memory = "none", era_items = "none")[[n]],
       cells = lapply(seq_along(subs), function(i) list(subject = subs[i], status = status_of(n)[i])))
})

save_json(list(
  example = TRUE,
  note = "Example project for designing the studio. Project files, explain and check are real rmvpadsl output; result maps and numbers are synthetic.",
  name = "objects-study", root = "~/studies/objects-study", version = proj$version,
  data = list(bids = "derivatives/lss", space = "MNI152NLin2009cAsym", mask = "brain"),
  sets = lapply(names(proj$config$data$sets), function(s) {
    st <- proj$config$data$sets[[s]]
    list(name = s, entities = st[setdiff(names(st), "keep")], keep = st$keep, subjects = length(subs), runs = 6L,
         rows_per_subject = 24L * 6L)
  }),
  folds = proj$config$folds, run = proj$config$run,
  subjects = subs,
  project_file = readLines(file.path(root, "rmvpa.yml")),
  analyses = analyses,
  check = chk
), "project.json")

vol <- function(x) { x[is.na(x)] <- NaN; b64(as.numeric(x)) }
save_json(list(dims = dims, spacing = c(4, 4, 4), template = jsonlite::base64_enc(as.raw(tpl8)),
               mask = jsonlite::base64_enc(as.raw(as.integer(mask)))), "template.json")
save_json(list(
  analysis = "decode_category", metric = "Accuracy", chance = chance, classes = 4L, metrics = c("Accuracy", "AUC"),
  null = list(kind = "searchlight_permutation", permutations = nperm, values = null_vals, centres_sampled = 0.1,
              n_null = n_null, p_floor_pooled = signif(1 / (n_null + 1), 3), p_resolution = round(1 / (nperm + 1), 4), correction = "fdr",
              level = "subject", adjust = "binned by searchlight size (number of voxels), so edge searchlights are compared with equally small ones",
              hist = list(edges = edges, counts = null_hist),
              description = "pooled searchlight permutation null, adjusted by searchlight size; not a voxelwise exact test"),
  centres = list(n = n_centres, values = round(sample(obs, min(6000, length(obs))), 4),
                 hist = list(edges = edges, counts = as.integer(table(cut(pmin(pmax(obs, 0.1), 0.4999), edges, include.lowest = TRUE))))),
  fdr = list(q = 0.05, threshold = round(fdr_thr, 4), surviving = sum(q < 0.05)),
  peak_z = kz - 1L,
  group = list(roi_rule = "centres where at least k of the other subjects survive their own FDR (q < .05)", roi_k = 3L,
               roi_sizes = as.list(setNames(round(rnorm(15, 160, 25)), subs[1:15])),
               test = NULL, note = "no group-level test is computed; these are per-subject values"),
  peak = list(ijk = as.integer(arrayInd(peak_idx, dims)) - 1L, value = round(grp[peak_idx], 4)),
  subjects = subj,
  maps = list(group_mean = vol(grp), prevalence = vol(prev), prevalence_all = vol(prev_all)),
  hist_edges = edges,
  sidecar = list(rmvpa_dsl = "0.1", rmvpa_version = "0.1.3", rmvpadsl_version = "0.1.0.9000", analysis = "decode_category",
                 subject = "01", verb = "decode", model_class = "mvpa_model",
                 inference = list(test = "permutation", permutations = nperm, shuffle = "within_run", null_method = "adjusted",
                                  adjust_by = "nfeatures", subsample = 0.1, correction = "fdr", p_floor = round(1 / (nperm + 1), 4)),
                 chance = list(Accuracy = chance, AUC = 0, classes = 4L), lock = "match",
                 rows = list(train = list(before_keep = 144L, kept = 131L)),
                 spec_hash = "9f2c41d0b7e3a5c6", folds = list(by = "run", scheme = "leave_one_out"))
), "results-decode_category.json")
save_json(list(analysis = "rsa_animacy", metric = "t", chance = 0, metrics = c("animacy", "category", "time"),
               null = list(kind = "searchlight_permutation", permutations = 1000L, joint = TRUE,
                           description = "joint null of no association with any predictor; per-predictor maps have no conditional test"),
               maps = list(animacy = vol(anim_t), category = vol(cat_t), time = vol(time_t))), "results-rsa_animacy.json")
save_json(list(analysis = "transfer_memory", metric = "Accuracy", chance = chance,
               chance_source = "1 / 4 classes; the retrieval set is balanced (36 trials per class per subject)",
               regions = xfer_tab, subjects = subs,
               atlas = list(file = "atlas/space-MNI152NLin2009cAsym_atlas-demo_dseg.nii.gz", md5 = unname(tools::md5sum(file.path(root, "atlas", "space-MNI152NLin2009cAsym_atlas-demo_dseg.nii.gz")))),
               null = list(kind = "none", description = "none: descriptive only (no permutation path in rMVPA for this analysis)")),
          "results-transfer_memory.json")
cat("Wrote", out, "\n")
