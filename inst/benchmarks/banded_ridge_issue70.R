# Reproducible simulation corresponding to GitHub issue #70.
#
# Run from the package root with, for example:
#   Rscript inst/benchmarks/banded_ridge_issue70.R 30
#
# The output reports absolute paired changes in outer-OOF correlation. It is a
# reproducibility/characterization artifact, not a superiority test.

issue70_ar1_matrix <- function(n, p, blocks, rho = 0.6) {
  out <- matrix(NA_real_, n, p)
  for (jj in seq_len(p)) {
    for (block in unique(blocks)) {
      rows <- which(blocks == block)
      innovations <- stats::rnorm(length(rows))
      values <- numeric(length(rows))
      values[[1L]] <- innovations[[1L]]
      if (length(rows) > 1L) {
        for (ii in 2:length(rows)) {
          values[[ii]] <- rho * values[[ii - 1L]] +
            sqrt(1 - rho^2) * innovations[[ii]]
        }
      }
      out[rows, jj] <- values
    }
  }
  out
}

issue70_dataset <- function(Y, dims) {
  n <- nrow(Y)
  data_array <- array(as.vector(t(Y)), c(dims, n))
  data_space <- neuroim2::NeuroSpace(c(dims, n), c(1, 1, 1))
  mask_space <- neuroim2::NeuroSpace(dims, c(1, 1, 1))
  rMVPA::mvpa_dataset(
    neuroim2::NeuroVec(data_array, data_space),
    mask = neuroim2::NeuroVol(array(1, dims), mask_space)
  )
}

issue70_one_replication <- function(seed,
                                    n = 120L,
                                    band_sizes = c(signal = 6L, nuisance = 4L),
                                    snr = rep(c(1, 0.4, 0), c(14L, 13L, 13L)),
                                    outer_folds = 6L,
                                    inner_folds = 5L,
                                    alphas = 10^seq(-2, 2, length.out = 5L),
                                    theta_grid_points = 4L) {
  stopifnot(length(snr) == 40L, sum(band_sizes) == 10L)
  set.seed(seed)
  if (n %% outer_folds != 0L) {
    stop("n must be divisible by outer_folds for contiguous segment blocks.",
         call. = FALSE)
  }
  blocks <- rep(seq_len(outer_folds), each = n %/% outer_folds)
  X <- issue70_ar1_matrix(n, sum(band_sizes), blocks)
  beta <- matrix(stats::rnorm(band_sizes[["signal"]] * length(snr)),
                 band_sizes[["signal"]], length(snr))
  beta <- sweep(beta, 2L, sqrt(colSums(beta^2)), "/")
  latent <- X[, seq_len(band_sizes[["signal"]]), drop = FALSE] %*% beta
  for (block in unique(blocks)) {
    rows <- which(blocks == block)
    latent[rows, ] <- sweep(latent[rows, , drop = FALSE], 2L,
                            colMeans(latent[rows, , drop = FALSE]), "-")
  }
  latent <- scale(latent)
  noise <- matrix(stats::rnorm(n * length(snr)), n, length(snr))
  for (block in unique(blocks)) {
    rows <- which(blocks == block)
    noise[rows, ] <- sweep(noise[rows, , drop = FALSE], 2L,
                           colMeans(noise[rows, , drop = FALSE]), "-")
  }
  Y <- sweep(latent, 2L, snr, "*") + noise

  features <- rMVPA::feature_sets(
    X, rMVPA::blocks(signal = band_sizes[["signal"]],
                     nuisance = band_sizes[["nuisance"]])
  )
  design <- rMVPA::feature_sets_design(
    features, block_var_train = blocks, time_series = TRUE
  )
  dataset <- issue70_dataset(Y, c(5L, 4L, 2L))

  make_model <- function(alpha_scope, theta_scope) {
    rMVPA::banded_ridge_model(
      dataset, design,
      outer_crossval = outer_folds,
      tune_crossval = inner_folds,
      alphas = alphas,
      theta_method = "grid",
      theta_grid_points = theta_grid_points,
      alpha_scope = alpha_scope,
      theta_scope = theta_scope,
      solver = "auto",
      target_batch_size = length(snr),
      seed = seed
    )
  }

  shared_model <- make_model("shared", "shared")
  response_model <- make_model("response", "response")
  shared_time <- system.time(shared <- rMVPA::run_banded_ridge(shared_model))
  response_time <- system.time(response <- rMVPA::run_banded_ridge(response_model))
  fold_manifest <- paste(vapply(shared_model$outer_folds, function(fold) {
    paste0(
      fold$id, "[test_rows=", min(fold$test), "-", max(fold$test),
      ";blocks=", paste(unique(blocks[fold$test]), collapse = ":"), "]"
    )
  }, character(1)), collapse = "|")
  candidate_manifest <- paste(shared_model$candidates$candidate_id,
                              collapse = ";")

  data.frame(
    replication = seq_along(snr) * 0L + seed,
    seed = seed,
    response = shared$metrics$response,
    voxel_index = shared$metrics$voxel_index,
    snr = snr,
    shared_oof_correlation = shared$metrics$correlation,
    response_oof_correlation = response$metrics$correlation,
    absolute_correlation_gain = response$metrics$correlation -
      shared$metrics$correlation,
    shared_oof_r2 = shared$metrics$r2,
    response_oof_r2 = response$metrics$r2,
    shared_elapsed_seconds = unname(shared_time[["elapsed"]]),
    response_elapsed_seconds = unname(response_time[["elapsed"]]),
    solver = shared$provenance$solver,
    fold_manifest = fold_manifest,
    candidate_manifest = candidate_manifest,
    stringsAsFactors = FALSE
  )
}

issue70_run <- function(seeds = 2026082801L + seq_len(30L),
                        output_file = file.path(
                          "inst", "extdata", "banded_ridge_issue70_results.csv"
                        )) {
  rows <- lapply(seeds, issue70_one_replication)
  out <- do.call(rbind, rows)
  out$package_version <- as.character(utils::packageVersion("rMVPA"))
  out$r_version <- R.version.string
  out$platform <- R.version$platform
  utils::write.csv(out, output_file, row.names = FALSE)
  out
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_replications <- if (length(args)) as.integer(args[[1L]]) else 30L
  output_file <- if (length(args) >= 2L) args[[2L]] else file.path(
    "inst", "extdata", "banded_ridge_issue70_results.csv"
  )
  result <- issue70_run(
    seeds = 2026082801L + seq_len(n_replications),
    output_file = output_file
  )
  summary_table <- aggregate(
    cbind(shared_oof_correlation, response_oof_correlation,
          absolute_correlation_gain) ~ snr,
    data = result, FUN = mean
  )
  print(summary_table, row.names = FALSE)
}
