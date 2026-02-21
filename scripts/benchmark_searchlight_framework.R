suppressPackageStartupMessages({
  library(devtools)
})

devtools::load_all(".", quiet = TRUE)

if (requireNamespace("futile.logger", quietly = TRUE)) {
  futile.logger::flog.threshold(futile.logger::ERROR)
}

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

flag_enabled <- function(name) {
  tolower(Sys.getenv(name, "false")) %in% c("1", "true", "yes", "on")
}

env_int <- function(name, default) {
  val <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (is.na(val) || val < 1L) default else val
}

env_csv <- function(name) {
  raw <- trimws(Sys.getenv(name, ""))
  if (!nzchar(raw)) {
    return(character())
  }
  vals <- trimws(unlist(strsplit(raw, ",", fixed = TRUE)))
  vals[nzchar(vals)]
}

resolve_bench_profile <- function(profile) {
  match.arg(tolower(profile), c("legacy", "phase1", "phase2"))
}

bench_profile_options <- function(profile) {
  profile <- resolve_bench_profile(profile)
  list(
    rMVPA.fold_cache_enabled = profile %in% c("phase1", "phase2"),
    rMVPA.matrix_first_roi = identical(profile, "phase2"),
    rMVPA.fast_filter_roi = identical(profile, "phase2"),
    rMVPA.naive_xdec_fast_kernel = identical(profile, "phase2")
  )
}

with_profile_timing <- function(expr) {
  old <- options(rMVPA.profile_searchlight = TRUE)
  on.exit(options(old), add = TRUE)
  force(expr)
}

build_mvpa_spec <- function() {
  ds <- gen_sample_dataset(D = c(7, 7, 7), nobs = 60, nlevels = 3, blocks = 5)
  cv <- blocked_cross_validation(ds$design$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
}

build_rsa_spec <- function() {
  ds <- gen_sample_dataset(D = c(7, 7, 7), nobs = 48, nlevels = 3, blocks = 4)
  dmat <- dist(matrix(rnorm(48 * 48), 48, 48))
  rdes <- rsa_design(
    ~ dmat,
    list(dmat = dmat, block = ds$design$block_var),
    block_var = "block"
  )
  rsa_model(ds$dataset, rdes, regtype = "pearson")
}

build_vector_rsa_spec <- function() {
  ds <- gen_sample_dataset(D = c(7, 7, 7), nobs = 48, nlevels = 4, blocks = 4)
  label_levels <- letters[1:8]
  labels <- factor(rep(label_levels, length.out = 48), levels = label_levels)
  dmat <- as.matrix(dist(matrix(rnorm(length(label_levels) * 12), length(label_levels), 12)))
  rownames(dmat) <- label_levels
  colnames(dmat) <- label_levels
  vdes <- vector_rsa_design(
    D = dmat,
    labels = labels,
    block_var = ds$design$block_var
  )
  vector_rsa_model(ds$dataset, vdes, rsa_simfun = "pearson")
}

build_naive_xdec_spec <- function() {
  ds <- gen_sample_dataset(
    D = c(7, 7, 7),
    nobs = 60,
    nlevels = 3,
    blocks = 5,
    external_test = TRUE,
    ntest_obs = 60
  )
  naive_xdec_model(ds$dataset, ds$design, return_predictions = FALSE)
}

build_clustered_spec <- function() {
  ds <- gen_clustered_sample_dataset(D = c(8, 8, 8), nobs = 40, K = 10, blocks = 4, nlevels = 3)
  cv <- blocked_cross_validation(ds$design$block_var)
  mvpa_model(
    model = load_model("sda_notune"),
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv
  )
}

collect_row <- function(workload, repetition, elapsed, res, backend, bench_profile) {
  timing <- attr(res, "timing")
  iter_timing <- timing$iterate %||% list()
  iter_totals <- iter_timing$totals %||% list()
  bad_results <- attr(res, "bad_results")
  n_bad <- if (is.data.frame(bad_results)) nrow(bad_results) else NA_integer_

  data.frame(
    workload = workload,
    repetition = repetition,
    backend = backend,
    bench_profile = bench_profile,
    fold_cache_enabled = isTRUE(getOption("rMVPA.fold_cache_enabled", FALSE)),
    matrix_first_roi = isTRUE(getOption("rMVPA.matrix_first_roi", FALSE)),
    fast_filter_roi = isTRUE(getOption("rMVPA.fast_filter_roi", FALSE)),
    naive_xdec_fast_kernel = isTRUE(getOption("rMVPA.naive_xdec_fast_kernel", FALSE)),
    elapsed_seconds = as.numeric(elapsed),
    setup_seconds = as.numeric(timing$setup_seconds %||% NA_real_),
    iterate_seconds = as.numeric(timing$iterate_seconds %||% NA_real_),
    combine_seconds = as.numeric(timing$combine_seconds %||% NA_real_),
    iter_get_samples_seconds = as.numeric(iter_totals$get_samples_seconds %||% NA_real_),
    iter_roi_extract_seconds = as.numeric(iter_totals$roi_extract_seconds %||% NA_real_),
    iter_run_future_seconds = as.numeric(iter_totals$run_future_seconds %||% NA_real_),
    iter_batch_seconds = as.numeric(iter_totals$batch_seconds %||% NA_real_),
    n_centers = as.integer(timing$n_centers %||% NA_integer_),
    n_good_results = as.integer(timing$n_good_results %||% NA_integer_),
    n_bad_results = as.integer(timing$n_bad_results %||% n_bad),
    stringsAsFactors = FALSE
  )
}

run_case <- function(workload, builder, radius, nrep, backend, bench_profile, ..., warmup = TRUE) {
  cat(sprintf("\n[%s] building model spec...\n", workload))
  mspec <- builder()

  if (isTRUE(warmup)) {
    invisible(run_searchlight(mspec, radius = radius, method = "standard", backend = backend, ...))
  }

  rows <- vector("list", nrep)
  for (i in seq_len(nrep)) {
    cat(sprintf("[%s] repetition %d/%d\n", workload, i, nrep))
    elapsed <- NA_real_
    res <- NULL

    out <- tryCatch({
      elapsed <- as.numeric(system.time({
        res <- with_profile_timing(
          run_searchlight(mspec, radius = radius, method = "standard", backend = backend, ...)
        )
      })["elapsed"])
      collect_row(workload, i, elapsed, res, backend, bench_profile = bench_profile)
    }, error = function(e) {
      data.frame(
        workload = workload,
        repetition = i,
        backend = backend,
        bench_profile = bench_profile,
        fold_cache_enabled = isTRUE(getOption("rMVPA.fold_cache_enabled", FALSE)),
        matrix_first_roi = isTRUE(getOption("rMVPA.matrix_first_roi", FALSE)),
        fast_filter_roi = isTRUE(getOption("rMVPA.fast_filter_roi", FALSE)),
        naive_xdec_fast_kernel = isTRUE(getOption("rMVPA.naive_xdec_fast_kernel", FALSE)),
        elapsed_seconds = NA_real_,
        setup_seconds = NA_real_,
        iterate_seconds = NA_real_,
        combine_seconds = NA_real_,
        iter_get_samples_seconds = NA_real_,
        iter_roi_extract_seconds = NA_real_,
        iter_run_future_seconds = NA_real_,
        iter_batch_seconds = NA_real_,
        n_centers = NA_integer_,
        n_good_results = NA_integer_,
        n_bad_results = NA_integer_,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    })
    rows[[i]] <- out
  }

  do.call(rbind, rows)
}

summarise_rows <- function(rows) {
  split_rows <- split(rows, rows$workload)
  out <- lapply(split_rows, function(df) {
    data.frame(
      workload = df$workload[1],
      backend = df$backend[1],
      nrep = nrow(df),
      median_elapsed_seconds = stats::median(df$elapsed_seconds, na.rm = TRUE),
      median_setup_seconds = stats::median(df$setup_seconds, na.rm = TRUE),
      median_iterate_seconds = stats::median(df$iterate_seconds, na.rm = TRUE),
      median_combine_seconds = stats::median(df$combine_seconds, na.rm = TRUE),
      median_iter_get_samples_seconds = stats::median(df$iter_get_samples_seconds, na.rm = TRUE),
      median_iter_roi_extract_seconds = stats::median(df$iter_roi_extract_seconds, na.rm = TRUE),
      median_iter_run_future_seconds = stats::median(df$iter_run_future_seconds, na.rm = TRUE),
      median_iter_batch_seconds = stats::median(df$iter_batch_seconds, na.rm = TRUE),
      median_n_centers = stats::median(df$n_centers, na.rm = TRUE),
      median_n_bad_results = stats::median(df$n_bad_results, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

set.seed(9101)
nrep <- env_int("RMVPA_BENCH_REP", 3L)
backend <- Sys.getenv("RMVPA_BENCH_BACKEND", "default")
bench_profile <- resolve_bench_profile(Sys.getenv("RMVPA_BENCH_PROFILE", "legacy"))
include_clustered <- flag_enabled("RMVPA_BENCH_INCLUDE_CLUSTERED")
selected_cases <- env_csv("RMVPA_BENCH_CASES")
dry_run <- flag_enabled("RMVPA_BENCH_DRY_RUN")
options(bench_profile_options(bench_profile))

cases <- list(
  list(name = "mvpa_classification", builder = build_mvpa_spec, radius = 3, dots = list()),
  list(name = "rsa_model", builder = build_rsa_spec, radius = 3, dots = list()),
  list(name = "vector_rsa_model", builder = build_vector_rsa_spec, radius = 3, dots = list()),
  list(name = "naive_xdec_model", builder = build_naive_xdec_spec, radius = 3, dots = list())
)

if (isTRUE(include_clustered)) {
  cases[[length(cases) + 1L]] <- list(
    name = "clustered_mvpa",
    builder = build_clustered_spec,
    radius = 2,
    dots = list(k = 6)
  )
}

if (length(selected_cases) > 0L) {
  cases <- Filter(function(case) case$name %in% selected_cases, cases)
}
if (length(cases) == 0L) {
  stop("No benchmark cases selected. Check RMVPA_BENCH_CASES values.", call. = FALSE)
}

if (isTRUE(dry_run)) {
  manifest <- data.frame(
    workload = vapply(cases, function(case) case$name, character(1)),
    backend = backend,
    bench_profile = bench_profile,
    fold_cache_enabled = isTRUE(getOption("rMVPA.fold_cache_enabled", FALSE)),
    matrix_first_roi = isTRUE(getOption("rMVPA.matrix_first_roi", FALSE)),
    fast_filter_roi = isTRUE(getOption("rMVPA.fast_filter_roi", FALSE)),
    naive_xdec_fast_kernel = isTRUE(getOption("rMVPA.naive_xdec_fast_kernel", FALSE)),
    nrep = nrep,
    radius = vapply(cases, function(case) case$radius, numeric(1)),
    stringsAsFactors = FALSE
  )

  cat("\n=== Searchlight Framework Benchmark (Dry Run) ===\n")
  print(manifest, row.names = FALSE)

  out_summary <- Sys.getenv("RMVPA_BENCH_OUT", "")
  if (nzchar(out_summary)) {
    utils::write.csv(manifest, out_summary, row.names = FALSE)
    cat(sprintf("\nWrote summary CSV: %s\n", out_summary))
  }

  out_raw <- Sys.getenv("RMVPA_BENCH_OUT_RAW", "")
  if (nzchar(out_raw)) {
    utils::write.csv(manifest, out_raw, row.names = FALSE)
    cat(sprintf("Wrote raw CSV: %s\n", out_raw))
  }

  quit(save = "no", status = 0L)
}

raw_rows <- do.call(rbind, lapply(cases, function(case) {
  do.call(
    run_case,
    c(
      list(
        workload = case$name,
        builder = case$builder,
        radius = case$radius,
        nrep = nrep,
        backend = backend,
        bench_profile = bench_profile
      ),
      case$dots
    )
  )
}))

summary_rows <- summarise_rows(raw_rows)

cat("\n=== Searchlight Framework Benchmark (Raw) ===\n")
print(raw_rows, row.names = FALSE)

cat("\n=== Searchlight Framework Benchmark (Summary) ===\n")
print(summary_rows, row.names = FALSE)

out_summary <- Sys.getenv("RMVPA_BENCH_OUT", "")
if (nzchar(out_summary)) {
  utils::write.csv(summary_rows, out_summary, row.names = FALSE)
  cat(sprintf("\nWrote summary CSV: %s\n", out_summary))
}

out_raw <- Sys.getenv("RMVPA_BENCH_OUT_RAW", "")
if (nzchar(out_raw)) {
  utils::write.csv(raw_rows, out_raw, row.names = FALSE)
  cat(sprintf("Wrote raw CSV: %s\n", out_raw))
}
