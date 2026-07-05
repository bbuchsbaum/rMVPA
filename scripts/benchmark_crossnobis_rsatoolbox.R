#!/usr/bin/env Rscript

# Benchmarks rMVPA crossnobis on matched deterministic data.
#
# The report separates two comparisons:
# 1. Public API: rMVPA end-to-end from matrix + CV spec vs rsatoolbox public Dataset API.
# 2. Kernel: rMVPA distance kernel vs a NumPy implementation of the rsatoolbox formula
#    on precomputed fold means. This formula row is a reference baseline, not public API.

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
python_helper <- file.path(repo_root, "scripts", "benchmark_crossnobis_rsatoolbox.py")

if (!dir.exists(vendor_src)) {
  stop("Vendored rsatoolbox source not found at `vendor/rsatoolbox/src`.", call. = FALSE)
}
if (!file.exists(python_helper)) {
  stop("Python helper not found: ", python_helper, call. = FALSE)
}

compute_cv <- get("compute_crossvalidated_means_sl", envir = asNamespace("rMVPA"))
compute_dist <- get("compute_crossnobis_distances_sl", envir = asNamespace("rMVPA"))
compute_gram <- get("compute_crossnobis_second_moment_sl", envir = asNamespace("rMVPA"))

parse_cases <- function(spec) {
  if (!nzchar(spec)) {
    return(NULL)
  }
  pieces <- strsplit(spec, ",", fixed = TRUE)[[1]]
  lapply(pieces, function(piece) {
    vals <- as.integer(strsplit(piece, "x", fixed = TRUE)[[1]])
    if (length(vals) != 4 || anyNA(vals) || any(vals <= 0)) {
      stop("Invalid case spec `", piece, "`. Expected KxVxMxreps, e.g. 24x128x6x2.",
           call. = FALSE)
    }
    list(K = vals[[1]], V = vals[[2]], M = vals[[3]], reps = vals[[4]])
  })
}

default_cases <- function(profile) {
  switch(profile,
    quick = list(
      list(K = 12L, V = 64L, M = 4L, reps = 2L),
      list(K = 24L, V = 128L, M = 6L, reps = 2L),
      list(K = 40L, V = 200L, M = 8L, reps = 2L)
    ),
    stress = list(
      list(K = 24L, V = 128L, M = 6L, reps = 3L),
      list(K = 48L, V = 256L, M = 8L, reps = 2L),
      list(K = 80L, V = 512L, M = 8L, reps = 2L)
    ),
    stop("Unknown RMVPA_RSATB_BENCH_PROFILE `", profile, "`. Use `quick` or `stress`.",
         call. = FALSE)
  )
}

make_case_id <- function(case) {
  paste0("K", case$K, "_V", case$V, "_M", case$M, "_R", case$reps)
}

make_bench_data <- function(K, V, M, reps) {
  n_obs <- K * M * reps
  sl_data <- matrix(NA_real_, nrow = n_obs, ncol = V)
  cond <- character(n_obs)
  folds <- integer(n_obs)
  vox <- seq_len(V)

  row <- 1L
  for (fold in seq_len(M)) {
    for (rep in seq_len(reps)) {
      for (k in seq_len(K)) {
        sl_data[row, ] <- sin(0.07 * k * vox + 0.11 * fold) +
          cos(0.03 * rep * vox + 0.17 * k) +
          0.02 * fold - 0.015 * k
        cond[[row]] <- sprintf("C%03d", k)
        folds[[row]] <- fold
        row <- row + 1L
      }
    }
  }

  colnames(sl_data) <- sprintf("V%03d", seq_len(V))
  levels <- sprintf("C%03d", seq_len(K))
  list(
    sl_data = sl_data,
    cond = factor(cond, levels = levels),
    folds = folds,
    levels = levels
  )
}

summarise_microbenchmark <- function(mb, case) {
  df <- as.data.frame(mb)
  split_ms <- split(df$time / 1e6, as.character(df$expr))
  out <- do.call(rbind, lapply(names(split_ms), function(method) {
    vals <- split_ms[[method]]
    data.frame(
      case_id = make_case_id(case),
      method = method,
      K = case$K,
      V = case$V,
      M = case$M,
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
  dat <- make_bench_data(case$K, case$V, case$M, case$reps)
  mvpa_des <- structure(
    list(
      Y = dat$cond,
      design_matrix = matrix(0, nrow = nrow(dat$sl_data), ncol = 1)
    ),
    class = c("mvpa_design", "list")
  )
  cv_spec <- blocked_cross_validation(dat$folds)

  cv_outputs <- compute_cv(
    dat$sl_data, mvpa_des, cv_spec,
    estimation_method = "crossnobis",
    whitening_matrix_W = NULL,
    return_folds = TRUE
  )
  U_folds <- cv_outputs$fold_estimates
  r_dist <- compute_dist(U_folds, P_voxels = case$V)

  mb <- microbenchmark::microbenchmark(
    rmvpa_kernel_distance = compute_dist(U_folds, P_voxels = case$V),
    rmvpa_kernel_second_moment = compute_gram(U_folds, vectorize = FALSE),
    rmvpa_end_to_end = {
      cv_outputs_i <- compute_cv(
        dat$sl_data, mvpa_des, cv_spec,
        estimation_method = "crossnobis",
        whitening_matrix_W = NULL,
        return_folds = TRUE
      )
      compute_dist(cv_outputs_i$fold_estimates, P_voxels = case$V)
    },
    times = times
  )

  list(
    summary = summarise_microbenchmark(mb, case),
    distances = unname(r_dist)
  )
}

profile <- Sys.getenv("RMVPA_RSATB_BENCH_PROFILE", "quick")
cases <- parse_cases(Sys.getenv("RMVPA_RSATB_BENCH_CASES", ""))
if (is.null(cases)) {
  cases <- default_cases(profile)
}
for (i in seq_along(cases)) {
  cases[[i]]$case_id <- make_case_id(cases[[i]])
}

times <- as.integer(Sys.getenv("RMVPA_RSATB_BENCH_TIMES", "10"))
if (is.na(times) || times < 1L) {
  stop("RMVPA_RSATB_BENCH_TIMES must be a positive integer.", call. = FALSE)
}
warmup <- as.integer(Sys.getenv("RMVPA_RSATB_BENCH_WARMUP", "2"))
if (is.na(warmup) || warmup < 0L) {
  stop("RMVPA_RSATB_BENCH_WARMUP must be a non-negative integer.", call. = FALSE)
}

cat("Benchmarking crossnobis against vendored rsatoolbox\n")
cat("Profile:", profile, " Iterations:", times, " Warmup:", warmup, "\n")
cat("Cases:", paste(vapply(cases, `[[`, character(1), "case_id"), collapse = ", "), "\n\n")

r_results <- vector("list", length(cases))
for (i in seq_along(cases)) {
  cat("Running rMVPA:", cases[[i]]$case_id, "\n")
  r_results[[i]] <- bench_rmvpa_case(cases[[i]], times)
}

config_file <- tempfile("rsatoolbox-bench-config-", fileext = ".json")
python_out <- tempfile("rsatoolbox-bench-output-", fileext = ".json")
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

parity <- do.call(rbind, lapply(seq_along(cases), function(i) {
  case_id <- cases[[i]]$case_id
  py_idx <- which(py_cases$case_id == case_id)
  if (length(py_idx) == 0L) {
    stop("Missing rsatoolbox result for case ", case_id, call. = FALSE)
  }
  do.call(rbind, lapply(py_idx, function(idx) {
    diff <- max(abs(r_results[[i]]$distances - unlist(py_cases$dissimilarities[[idx]])))
    data.frame(
      case_id = case_id,
      method = py_cases$method[[idx]],
      max_abs_diff = diff,
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

rsat_median <- setNames(
  py_summary$median_ms[py_summary$method == "rsatoolbox_public"],
  py_summary$case_id[py_summary$method == "rsatoolbox_public"]
)
formula_median <- setNames(
  py_summary$median_ms[py_summary$method == "rsatoolbox_numpy_formula"],
  py_summary$case_id[py_summary$method == "rsatoolbox_numpy_formula"]
)

public_methods <- all_summary$method %in% c("rmvpa_end_to_end", "rsatoolbox_public")
kernel_methods <- all_summary$method %in% c("rmvpa_kernel_distance", "rsatoolbox_numpy_formula")

all_summary$public_api_speedup_vs_rsatoolbox <- NA_real_
all_summary$public_api_speedup_vs_rsatoolbox[public_methods] <-
  rsat_median[all_summary$case_id[public_methods]] / all_summary$median_ms[public_methods]

all_summary$kernel_speedup_vs_formula_reference <- NA_real_
all_summary$kernel_speedup_vs_formula_reference[kernel_methods] <-
  formula_median[all_summary$case_id[kernel_methods]] / all_summary$median_ms[kernel_methods]

all_summary <- all_summary[order(all_summary$case_id, all_summary$method), ]
rownames(all_summary) <- NULL

public_summary <- all_summary[all_summary$method %in% c("rmvpa_end_to_end", "rsatoolbox_public"), ]
kernel_summary <- all_summary[all_summary$method %in% c("rmvpa_kernel_distance", "rsatoolbox_numpy_formula"), ]
internal_summary <- all_summary[all_summary$method %in% c("rmvpa_kernel_second_moment"), ]

cat("\nParity max abs diff:\n")
print(parity, row.names = FALSE)

cat("\nComparison scope:\n")
cat("- Public API: `rmvpa_end_to_end` starts from a numeric matrix and CV spec; `rsatoolbox_public` starts from a Dataset and uses calc_rdm_crossnobis().\n")
cat("- Kernel: `rmvpa_kernel_distance` and `rsatoolbox_numpy_formula` both start from precomputed fold means. The NumPy formula is a reference implementation of the rsatoolbox algebra, not an official public API path.\n")
cat("- The vendored rsatoolbox Cython engine is stubbed only to allow imports; the balanced calc_rdm_crossnobis() path benchmarked here is pure Python.\n")

public_cols <- c("case_id", "method", "K", "V", "M", "reps", "iterations",
                 "median_ms", "min_ms", "max_ms", "mean_ms",
                 "public_api_speedup_vs_rsatoolbox")
kernel_cols <- c("case_id", "method", "K", "V", "M", "reps", "iterations",
                 "median_ms", "min_ms", "max_ms", "mean_ms",
                 "kernel_speedup_vs_formula_reference")
internal_cols <- c("case_id", "method", "K", "V", "M", "reps", "iterations",
                   "median_ms", "min_ms", "max_ms", "mean_ms")

cat("\nPublic API timing (milliseconds; speedup > 1 means faster than rsatoolbox public API):\n")
print(public_summary[, public_cols, drop = FALSE], row.names = FALSE, digits = 4)

cat("\nKernel timing (milliseconds; speedup > 1 means faster than the rsatoolbox formula reference):\n")
print(kernel_summary[, kernel_cols, drop = FALSE], row.names = FALSE, digits = 4)

cat("\nAdditional rMVPA internal kernel timing (not a direct rsatoolbox output comparison):\n")
print(internal_summary[, internal_cols, drop = FALSE], row.names = FALSE, digits = 4)

out_csv <- Sys.getenv("RMVPA_RSATB_BENCH_OUT", "")
if (nzchar(out_csv)) {
  utils::write.csv(all_summary, out_csv, row.names = FALSE)
  cat("\nWrote benchmark summary to", out_csv, "\n")
}
