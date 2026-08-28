# Fixed-shape public-path performance characterization for banded ridge.
#
# Run from the package root with:
#   Rscript inst/benchmarks/banded_ridge_performance.R
#
# Wall times are descriptive snapshots, not CI pass/fail thresholds. Structural
# allocation limits are enforced separately by the package and its tests.

brperf_hardware <- function() {
  chip <- NA_character_
  memory <- NA_character_
  if (identical(Sys.info()[["sysname"]], "Darwin") &&
      nzchar(Sys.which("system_profiler"))) {
    lines <- tryCatch(
      system2("system_profiler", "SPHardwareDataType", stdout = TRUE,
              stderr = FALSE),
      error = function(e) character()
    )
    chip_line <- grep("^[[:space:]]+Chip:", lines, value = TRUE)
    memory_line <- grep("^[[:space:]]+Memory:", lines, value = TRUE)
    if (length(chip_line)) chip <- sub(".*Chip:[[:space:]]*", "", chip_line[[1L]])
    if (length(memory_line)) {
      memory <- sub(".*Memory:[[:space:]]*", "", memory_line[[1L]])
    }
  }
  c(chip = chip, memory = memory)
}

brperf_dataset <- function(Y, dims) {
  n <- nrow(Y)
  data_array <- array(as.vector(t(Y)), c(dims, n))
  rMVPA::mvpa_dataset(
    neuroim2::NeuroVec(
      data_array, neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))
    ),
    mask = neuroim2::NeuroVol(
      array(1, dims), neuroim2::NeuroSpace(dims, c(1, 1, 1))
    )
  )
}

brperf_fixture <- function(seed, n, p, v = 12L) {
  stopifnot(v == 12L, p %% 2L == 0L, n %% 3L == 0L)
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  beta <- matrix(stats::rnorm(min(8L, p) * v), min(8L, p), v)
  Y <- X[, seq_len(nrow(beta)), drop = FALSE] %*% beta +
    matrix(stats::rnorm(n * v, sd = 0.5), n, v)
  features <- rMVPA::feature_sets(
    X, rMVPA::blocks(a = p %/% 2L, b = p %/% 2L)
  )
  design <- rMVPA::feature_sets_design(
    features, block_var_train = rep(1:3, each = n %/% 3L),
    time_series = TRUE
  )
  list(
    dataset = brperf_dataset(Y, c(3L, 2L, 2L)),
    design = design
  )
}

brperf_one_shape <- function(shape_id,
                             seed,
                             n,
                             p,
                             v = 12L,
                             iterations = 5L) {
  fixture <- brperf_fixture(seed, n, p, v)
  solvers <- c("direct", "svd_primal", "dual_kernel")
  outputs <- list()
  timings <- list()
  for (solver in solvers) {
    model <- rMVPA::banded_ridge_model(
      fixture$dataset, fixture$design,
      outer_crossval = 3L, tune_crossval = 2L,
      alphas = c(0.1, 1, 10), theta_method = "fixed",
      theta = c(a = 0.5, b = 0.5), solver = solver,
      target_batch_size = v, return_predictions = TRUE,
      seed = seed + 1L
    )
    elapsed <- numeric(iterations)
    for (ii in seq_len(iterations)) {
      timed <- system.time(outputs[[solver]] <- rMVPA::run_banded_ridge(model))
      elapsed[[ii]] <- unname(timed[["elapsed"]])
    }
    timings[[solver]] <- elapsed
  }
  reference <- outputs$direct$predictions
  do.call(rbind, lapply(solvers, function(solver) {
    out <- outputs[[solver]]
    plan <- out$provenance$solver_plan
    data.frame(
      shape = shape_id,
      seed = seed,
      n = n,
      p = p,
      v = v,
      outer_folds = 3L,
      inner_folds = 2L,
      theta_candidates = 1L,
      alpha_candidates = 3L,
      response_batch_size = v,
      solver = solver,
      elapsed_median_seconds = stats::median(timings[[solver]]),
      elapsed_min_seconds = min(timings[[solver]]),
      elapsed_max_seconds = max(timings[[solver]]),
      iterations = iterations,
      max_abs_prediction_error_vs_direct =
        max(abs(out$predictions - reference)),
      estimated_solver_intermediate_peak_bytes =
        unname(plan$estimated_bytes[[solver]]),
      max_chunk_result_bytes = out$provenance$max_chunk_result_bytes,
      retained_output_estimated_bytes =
        unname(out$provenance$retention_estimated_bytes[["total"]]),
      decomposition_count = out$provenance$solver_decomposition_count,
      cache_hits = out$provenance$solver_cache_hits,
      stringsAsFactors = FALSE
    )
  }))
}

brperf_run <- function(output_file = file.path(
                         "inst", "extdata",
                         "banded_ridge_performance_results.csv"
                       ),
                       iterations = 5L) {
  shapes <- list(
    small = list(seed = 2026082851L, n = 96L, p = 24L, v = 12L),
    wide = list(seed = 2026082852L, n = 60L, p = 512L, v = 12L)
  )
  rows <- do.call(rbind, lapply(names(shapes), function(id) {
    args <- shapes[[id]]
    do.call(brperf_one_shape, c(list(shape_id = id), args,
                                list(iterations = iterations)))
  }))
  hardware <- brperf_hardware()
  rows$package_version <- as.character(utils::packageVersion("rMVPA"))
  rows$r_version <- R.version.string
  rows$platform <- R.version$platform
  rows$hardware_chip <- unname(hardware[["chip"]])
  rows$hardware_memory <- unname(hardware[["memory"]])
  rows$logical_cores_reported <- parallel::detectCores(logical = TRUE)
  rows$git_head <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[[1L]],
    error = function(e) NA_character_
  )
  rows$worktree_dirty <- tryCatch(
    length(system2("git", c("status", "--porcelain"), stdout = TRUE,
                   stderr = FALSE)) > 0L,
    error = function(e) NA
  )
  utils::write.csv(rows, output_file, row.names = FALSE)
  rows
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_file <- if (length(args)) args[[1L]] else file.path(
    "inst", "extdata", "banded_ridge_performance_results.csv"
  )
  result <- brperf_run(output_file)
  print(result[, c(
    "shape", "solver", "elapsed_median_seconds",
    "estimated_solver_intermediate_peak_bytes",
    "max_abs_prediction_error_vs_direct"
  )], row.names = FALSE)
}
