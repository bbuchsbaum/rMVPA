suppressPackageStartupMessages({
  library(devtools)
})

devtools::load_all(".", quiet = TRUE)

set.seed(8801)

flag_enabled <- function(name) {
  tolower(Sys.getenv(name, "false")) %in% c("1", "true", "yes", "on")
}

build_spec <- function(D = c(7, 7, 7), nobs = 60, nlevels = 3, blocks = 5, gamma = 1e-2) {
  ds <- gen_sample_dataset(D = D, nobs = nobs, nlevels = nlevels, blocks = blocks)
  cv <- blocked_cross_validation(ds$design$block_var)
  mdl <- load_model("dual_lda")
  mvpa_model(
    model = mdl,
    dataset = ds$dataset,
    design = ds$design,
    model_type = "classification",
    crossval = cv,
    tune_grid = data.frame(gamma = gamma)
  )
}

bench_dual_lda <- function(mspec, radius = 3, nrep = 3L) {
  invisible(run_searchlight(mspec, radius = radius, method = "standard", incremental = TRUE))
  invisible(run_searchlight(mspec, radius = radius, method = "standard", incremental = FALSE))

  elapsed <- function(incremental) {
    as.numeric(system.time(
      run_searchlight(mspec, radius = radius, method = "standard", incremental = incremental)
    )["elapsed"])
  }

  inc_times <- replicate(nrep, elapsed(TRUE))
  full_times <- replicate(nrep, elapsed(FALSE))

  med_inc <- stats::median(inc_times)
  med_full <- stats::median(full_times)
  speedup <- med_full / med_inc

  data.frame(
    mode = c("incremental", "full_recompute"),
    median_seconds = c(med_inc, med_full),
    mean_seconds = c(mean(inc_times), mean(full_times)),
    stringsAsFactors = FALSE
  ) -> tbl

  cat("\nDual-LDA Searchlight Benchmark\n")
  cat(sprintf("Radius: %s | Repetitions: %d\n", radius, nrep))
  print(tbl, row.names = FALSE)
  cat(sprintf("\nEstimated speedup (full / incremental): %.2fx\n", speedup))
  invisible(list(summary = tbl, speedup = speedup, inc_times = inc_times, full_times = full_times))
}

spec <- build_spec()
bench_dual_lda(spec, radius = 3, nrep = 3L)

if (flag_enabled("RMVPA_BENCH_NIGHTLY")) {
  cat("\nRunning nightly large benchmark profile...\n")
  spec_nightly <- build_spec(D = c(10, 10, 10), nobs = 120, nlevels = 3, blocks = 6, gamma = 1e-2)
  bench_dual_lda(spec_nightly, radius = 4, nrep = 2L)
}
