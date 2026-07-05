#!/usr/bin/env Rscript

# Benchmarks ordinary RDM calculators against vendored rsatoolbox.
#
# The report separates two comparisons:
# 1. Public calculator path: rMVPA condition averaging + pairwise_dist() vs
#    rsatoolbox calc_rdm(..., descriptor=) for correlation and Euclidean RDMs.
# 2. Formula path: matched R and NumPy formulas on precomputed condition means,
#    including Mahalanobis with an explicit supplied precision matrix.

thread_env <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS")
for (var in thread_env) {
  if (!nzchar(Sys.getenv(var))) {
    do.call(Sys.setenv, as.list(stats::setNames("1", var)))
  }
}

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package `jsonlite` is required for this benchmark.", call. = FALSE)
}
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("Package `microbenchmark` is required for this benchmark.", call. = FALSE)
}

if (requireNamespace("pkgload", quietly = TRUE) && file.exists("DESCRIPTION")) {
  pkgload::load_all(".", quiet = TRUE)
} else {
  library(rMVPA)
}

repo_root <- normalizePath(".", mustWork = TRUE)
vendor_src <- file.path(repo_root, "vendor", "rsatoolbox", "src")
python_helper <- file.path(repo_root, "scripts", "benchmark_rdm_calculators_rsatoolbox.py")

if (!dir.exists(vendor_src)) {
  stop("Vendored rsatoolbox source not found at `vendor/rsatoolbox/src`.", call. = FALSE)
}
if (!file.exists(python_helper)) {
  stop("Python helper not found: ", python_helper, call. = FALSE)
}

pairwise_dist <- get("pairwise_dist", envir = asNamespace("rMVPA"))
group_means <- get("group_means", envir = asNamespace("rMVPA"))
cordist <- get("cordist", envir = asNamespace("rMVPA"))
eucdist <- get("eucdist", envir = asNamespace("rMVPA"))
rdm_vector_correlation <- get(".rdm_vector_correlation", envir = asNamespace("rMVPA"))
rdm_vector_sqeuclidean <- get(".rdm_vector_sqeuclidean", envir = asNamespace("rMVPA"))
rdm_vector_sqmahalanobis <- get(".rdm_vector_sqmahalanobis", envir = asNamespace("rMVPA"))

parse_cases <- function(spec) {
  if (!nzchar(spec)) {
    return(NULL)
  }
  pieces <- strsplit(spec, ",", fixed = TRUE)[[1]]
  lapply(pieces, function(piece) {
    vals <- as.integer(strsplit(piece, "x", fixed = TRUE)[[1]])
    if (length(vals) != 3 || anyNA(vals) || any(vals <= 0)) {
      stop("Invalid case spec `", piece, "`. Expected KxVxreps, e.g. 24x128x4.",
           call. = FALSE)
    }
    list(K = vals[[1]], V = vals[[2]], reps = vals[[3]])
  })
}

default_cases <- function(profile) {
  switch(profile,
    quick = list(
      list(K = 12L, V = 64L, reps = 3L),
      list(K = 24L, V = 128L, reps = 4L),
      list(K = 40L, V = 200L, reps = 4L)
    ),
    stress = list(
      list(K = 24L, V = 128L, reps = 6L),
      list(K = 48L, V = 256L, reps = 5L),
      list(K = 80L, V = 512L, reps = 4L)
    ),
    stop("Unknown RMVPA_RDM_BENCH_PROFILE `", profile, "`. Use `quick` or `stress`.",
         call. = FALSE)
  )
}

make_case_id <- function(case) {
  paste0("K", case$K, "_V", case$V, "_R", case$reps)
}

make_bench_data <- function(K, V, reps) {
  n_obs <- K * reps
  measurements <- matrix(NA_real_, nrow = n_obs, ncol = V)
  cond <- character(n_obs)
  vox <- seq_len(V)

  row <- 1L
  for (rep in seq_len(reps)) {
    for (k in seq_len(K)) {
      measurements[row, ] <- sin(0.07 * k * vox + 0.11 * rep) +
        cos(0.03 * rep * vox + 0.17 * k) +
        0.02 * rep -
        0.015 * k +
        0.001 * (k %% 3L) * vox / V
      cond[[row]] <- sprintf("C%03d", k)
      row <- row + 1L
    }
  }

  colnames(measurements) <- sprintf("V%03d", seq_len(V))
  levels <- sprintf("C%03d", seq_len(K))
  list(
    measurements = measurements,
    cond = factor(cond, levels = levels),
    levels = levels
  )
}

make_precision <- function(V) {
  idx <- seq_len(V)
  diag_vals <- 1 + (idx %% 7L) / 10
  u <- sin(0.17 * idx)
  diag(diag_vals, nrow = V, ncol = V) + 0.02 * tcrossprod(u) / V
}

condition_means <- function(measurements, cond, levels) {
  means <- group_means(measurements, margin = 1, group = factor(cond, levels = levels))
  means <- means[levels, , drop = FALSE]
  rownames(means) <- levels
  means
}

rdm_vec <- function(dmat) {
  unname(dmat[lower.tri(dmat)])
}

rdm_correlation_formula <- function(means) {
  rdm_vector_correlation(means, method = "pearson")
}

rdm_euclidean_formula <- function(means) {
  rdm_vector_sqeuclidean(means, normalize_by_features = TRUE)
}

rdm_mahalanobis_formula <- function(means, precision) {
  rdm_vector_sqmahalanobis(means, precision, normalize_by_features = TRUE)
}

rmvpa_pairwise_correlation <- function(measurements, cond, levels) {
  means <- condition_means(measurements, cond, levels)
  rdm_vec(pairwise_dist(cordist(method = "pearson"), means))
}

rmvpa_pairwise_euclidean <- function(measurements, cond, levels) {
  means <- condition_means(measurements, cond, levels)
  rdm_vec(pairwise_dist(eucdist(), means)^2) / ncol(means)
}

summarise_microbenchmark <- function(mb, case) {
  df <- as.data.frame(mb)
  split_ms <- split(df$time / 1e6, as.character(df$expr))
  out <- do.call(rbind, lapply(names(split_ms), function(expr) {
    vals <- split_ms[[expr]]
    distance <- sub("^.*_(correlation|euclidean|mahalanobis)$", "\\1", expr)
    method <- sub("_(correlation|euclidean|mahalanobis)$", "", expr)
    data.frame(
      case_id = make_case_id(case),
      method = method,
      distance = distance,
      K = case$K,
      V = case$V,
      reps = case$reps,
      iterations = length(vals),
      median_ms = median(vals),
      min_ms = min(vals),
      max_ms = max(vals),
      mean_ms = mean(vals),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

bench_rmvpa_case <- function(case, times) {
  dat <- make_bench_data(case$K, case$V, case$reps)
  means <- condition_means(dat$measurements, dat$cond, dat$levels)
  precision <- make_precision(case$V)

  r_distances <- list(
    correlation = list(
      rmvpa_pairwise_dist = rmvpa_pairwise_correlation(dat$measurements, dat$cond, dat$levels),
      rmvpa_r_formula = rdm_correlation_formula(means)
    ),
    euclidean = list(
      rmvpa_pairwise_dist = rmvpa_pairwise_euclidean(dat$measurements, dat$cond, dat$levels),
      rmvpa_r_formula = rdm_euclidean_formula(means)
    ),
    mahalanobis = list(
      rmvpa_r_formula = rdm_mahalanobis_formula(means, precision)
    )
  )

  mb <- microbenchmark::microbenchmark(
    rmvpa_pairwise_dist_correlation = rmvpa_pairwise_correlation(dat$measurements, dat$cond, dat$levels),
    rmvpa_pairwise_dist_euclidean = rmvpa_pairwise_euclidean(dat$measurements, dat$cond, dat$levels),
    rmvpa_r_formula_correlation = rdm_correlation_formula(means),
    rmvpa_r_formula_euclidean = rdm_euclidean_formula(means),
    rmvpa_r_formula_mahalanobis = rdm_mahalanobis_formula(means, precision),
    times = times
  )

  list(
    summary = summarise_microbenchmark(mb, case),
    distances = r_distances
  )
}

extract_py_distances <- function(py_cases, idx) {
  dissimilarities <- py_cases$dissimilarities
  if (is.matrix(dissimilarities)) {
    return(as.numeric(dissimilarities[idx, ]))
  }
  as.numeric(unlist(dissimilarities[[idx]], use.names = FALSE))
}

get_py_distances <- function(py_cases, case_id, method, distance) {
  idx <- which(py_cases$case_id == case_id &
                 py_cases$method == method &
                 py_cases$distance == distance)
  if (length(idx) != 1L) {
    stop("Expected one Python result for case=", case_id,
         " method=", method, " distance=", distance,
         "; got ", length(idx), ".", call. = FALSE)
  }
  extract_py_distances(py_cases, idx)
}

profile <- Sys.getenv("RMVPA_RDM_BENCH_PROFILE", "quick")
cases <- parse_cases(Sys.getenv("RMVPA_RDM_BENCH_CASES", ""))
if (is.null(cases)) {
  cases <- default_cases(profile)
}
for (i in seq_along(cases)) {
  cases[[i]]$case_id <- make_case_id(cases[[i]])
}

times <- as.integer(Sys.getenv("RMVPA_RDM_BENCH_TIMES", "10"))
if (is.na(times) || times < 1L) {
  stop("RMVPA_RDM_BENCH_TIMES must be a positive integer.", call. = FALSE)
}
warmup <- as.integer(Sys.getenv("RMVPA_RDM_BENCH_WARMUP", "2"))
if (is.na(warmup) || warmup < 0L) {
  stop("RMVPA_RDM_BENCH_WARMUP must be a non-negative integer.", call. = FALSE)
}

cat("Benchmarking ordinary RDM calculators against vendored rsatoolbox\n")
cat("Profile:", profile, " Iterations:", times, " Warmup:", warmup, "\n")
cat("Cases:", paste(vapply(cases, `[[`, character(1), "case_id"), collapse = ", "), "\n\n")

r_results <- vector("list", length(cases))
for (i in seq_along(cases)) {
  cat("Running rMVPA:", cases[[i]]$case_id, "\n")
  r_results[[i]] <- bench_rmvpa_case(cases[[i]], times)
}

config_file <- tempfile("rsatoolbox-rdm-bench-config-", fileext = ".json")
python_out <- tempfile("rsatoolbox-rdm-bench-output-", fileext = ".json")
jsonlite::write_json(
  list(vendor_src = vendor_src, times = times, warmup = warmup, cases = cases),
  config_file,
  auto_unbox = TRUE
)

python_bin <- Sys.getenv("RMVPA_PYTHON", "python")
cat("Running rsatoolbox via", python_bin, "\n")
status <- system2(
  python_bin,
  c(python_helper, "--config", config_file, "--output", python_out),
  stdout = TRUE,
  stderr = TRUE
)
exit_status <- attr(status, "status")
if (!is.null(exit_status) && exit_status != 0L) {
  cat(paste(status, collapse = "\n"), "\n")
  stop("rsatoolbox benchmark helper failed with status ", exit_status, call. = FALSE)
}
if (length(status) > 0L) {
  cat(paste(status, collapse = "\n"), "\n")
}

py_payload <- jsonlite::read_json(python_out, simplifyVector = TRUE)
py_cases <- py_payload$cases
py_summary <- py_cases[, setdiff(names(py_cases), "dissimilarities"), drop = FALSE]

parity_specs <- data.frame(
  distance = c("correlation", "correlation", "euclidean", "euclidean",
               "mahalanobis", "mahalanobis"),
  rmvpa_method = c("rmvpa_pairwise_dist", "rmvpa_r_formula",
                   "rmvpa_pairwise_dist", "rmvpa_r_formula",
                   "rmvpa_r_formula", "rmvpa_r_formula"),
  rsatoolbox_method = c("rsatoolbox_public", "rsatoolbox_numpy_formula",
                        "rsatoolbox_public", "rsatoolbox_numpy_formula",
                        "rsatoolbox_public", "rsatoolbox_numpy_formula"),
  comparison = c("public", "formula", "public", "formula",
                 "supplied_precision_public", "formula"),
  stringsAsFactors = FALSE
)

parity <- do.call(rbind, lapply(seq_along(cases), function(i) {
  case_id <- cases[[i]]$case_id
  do.call(rbind, lapply(seq_len(nrow(parity_specs)), function(j) {
    spec <- parity_specs[j, ]
    r_vec <- r_results[[i]]$distances[[spec$distance]][[spec$rmvpa_method]]
    py_vec <- get_py_distances(py_cases, case_id, spec$rsatoolbox_method, spec$distance)
    data.frame(
      case_id = case_id,
      distance = spec$distance,
      comparison = spec$comparison,
      rmvpa_method = spec$rmvpa_method,
      rsatoolbox_method = spec$rsatoolbox_method,
      max_abs_diff = max(abs(r_vec - py_vec)),
      stringsAsFactors = FALSE
    )
  }))
}))

if (any(parity$max_abs_diff > 1e-9)) {
  print(parity)
  stop("Parity check failed before speed comparison.", call. = FALSE)
}

r_summary <- do.call(rbind, lapply(r_results, `[[`, "summary"))
all_summary <- rbind(
  r_summary,
  py_summary[, names(r_summary), drop = FALSE]
)

all_summary$public_api_speedup_vs_rsatoolbox <- NA_real_
all_summary$formula_speedup_vs_numpy_reference <- NA_real_

for (case_id in unique(all_summary$case_id)) {
  for (distance in c("correlation", "euclidean")) {
    public_ref <- all_summary$median_ms[
      all_summary$case_id == case_id &
        all_summary$distance == distance &
        all_summary$method == "rsatoolbox_public"
    ]
    public_rows <- all_summary$case_id == case_id &
      all_summary$distance == distance &
      all_summary$method %in% c("rmvpa_pairwise_dist", "rsatoolbox_public")
    all_summary$public_api_speedup_vs_rsatoolbox[public_rows] <-
      public_ref / all_summary$median_ms[public_rows]
  }

  for (distance in c("correlation", "euclidean", "mahalanobis")) {
    formula_ref <- all_summary$median_ms[
      all_summary$case_id == case_id &
        all_summary$distance == distance &
        all_summary$method == "rsatoolbox_numpy_formula"
    ]
    formula_rows <- all_summary$case_id == case_id &
      all_summary$distance == distance &
      all_summary$method %in% c("rmvpa_r_formula", "rsatoolbox_numpy_formula")
    all_summary$formula_speedup_vs_numpy_reference[formula_rows] <-
      formula_ref / all_summary$median_ms[formula_rows]
  }
}

all_summary <- all_summary[order(all_summary$case_id, all_summary$distance, all_summary$method), ]
rownames(all_summary) <- NULL

public_summary <- all_summary[
  all_summary$distance %in% c("correlation", "euclidean") &
    all_summary$method %in% c("rmvpa_pairwise_dist", "rsatoolbox_public"),
]
formula_summary <- all_summary[
  all_summary$method %in% c("rmvpa_r_formula", "rsatoolbox_numpy_formula"),
]
mahalanobis_public_reference <- all_summary[
  all_summary$distance == "mahalanobis" &
    all_summary$method == "rsatoolbox_public",
]

cat("\nParity max abs diff:\n")
print(parity, row.names = FALSE)

cat("\nComparison scope:\n")
cat("- Public calculator path: `rmvpa_pairwise_dist` averages raw observations with rMVPA `group_means()` and then uses `pairwise_dist()`; `rsatoolbox_public` starts from a Dataset and uses calc_rdm(..., descriptor='cond').\n")
cat("- Euclidean parity uses rsatoolbox's convention: squared Euclidean distance divided by feature count.\n")
cat("- Mahalanobis uses an explicit supplied precision matrix and is compared as a formula/reference path because rMVPA's current `mahadist()` estimates shrinkage precision internally.\n")
cat("- The vendored rsatoolbox Cython engine is stubbed only to allow imports; these balanced RDM calculator paths do not depend on the unbuilt engine.\n")

public_cols <- c("case_id", "distance", "method", "K", "V", "reps", "iterations",
                 "median_ms", "min_ms", "max_ms", "mean_ms",
                 "public_api_speedup_vs_rsatoolbox")
formula_cols <- c("case_id", "distance", "method", "K", "V", "reps", "iterations",
                  "median_ms", "min_ms", "max_ms", "mean_ms",
                  "formula_speedup_vs_numpy_reference")
mahal_cols <- c("case_id", "distance", "method", "K", "V", "reps", "iterations",
                "median_ms", "min_ms", "max_ms", "mean_ms")

cat("\nPublic calculator timing (milliseconds; speedup > 1 means faster than rsatoolbox public API):\n")
print(public_summary[, public_cols, drop = FALSE], row.names = FALSE, digits = 4)

cat("\nFormula timing (milliseconds; speedup > 1 means faster than the NumPy formula reference):\n")
print(formula_summary[, formula_cols, drop = FALSE], row.names = FALSE, digits = 4)

cat("\nrsatoolbox public Mahalanobis timing with supplied precision (reference only; no current rMVPA public precision path):\n")
print(mahalanobis_public_reference[, mahal_cols, drop = FALSE], row.names = FALSE, digits = 4)

out_csv <- Sys.getenv("RMVPA_RDM_BENCH_OUT", "")
if (nzchar(out_csv)) {
  utils::write.csv(all_summary, out_csv, row.names = FALSE)
  cat("\nWrote benchmark summary to", out_csv, "\n")
}
